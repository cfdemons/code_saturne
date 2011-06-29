# -*- coding: utf-8 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2011 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module manages the differents possible outputs :
- listing printing
- post-processing and relationship with the FVM library
- monitoring points
- writer
- mesh

This module defines the following classes:
- OutputControlModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import string, sys, unittest
from types import FloatType

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

import Base.Toolbox as Tool
from Base.XMLvariables import Model
from Base.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class OutputControlModel(Model):
    """
    Class for Variables and Scalar model initialization.
    """
    def __init__(self, case):
        """
        Constructor
        """
        self.case = case
        node_control  = self.case.xmlGetNode('analysis_control')
        self.node_out = node_control.xmlInitNode('output')


    def defaultInitialValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['listing_printing_frequency'] = 1
        default['probe_recording_frequency'] = 1
        default['probe_recording_frequency_time'] = 0.1
        if self.case['salome']:
            default['postprocessing_format'] = "MED"
        default['probe_format'] = "DAT"
        default['coordinate'] = 0.0

        return default


    def __getCoordinates(self, name, coord):
        """
        Private method: return coordinate 'coord' for probe named 'name'
        """
        val = self.node_out.xmlGetNode('probe', name = name).xmlGetDouble(coord)
        if val == None:
            val = self.defaultInitialValues()['coordinate']
            self.__setCoordinates(name, coord, val)
        return val


    def __setCoordinates(self, name, coord, val):
        """
        Private method: put value of coordinate 'coord' for probe named 'name'
        """
        self.node_out.xmlGetNode('probe', name=name).xmlSetData(coord, val)


    def getListingFrequency(self):
        """
        Return the frequency for printing listing
        """
        f = self.node_out.xmlGetInt('listing_printing_frequency')
        if f == None:
            f = self.defaultInitialValues()['listing_printing_frequency']
            self.setListingFrequency(f)
        return f


    def setListingFrequency(self, freq):
        """
        Set the frequency for printing listing
        """
        self.isInt(freq)
        self.node_out.xmlSetData('listing_printing_frequency', freq)


    def defaultWriterValues(self):
        """Return the default values - Method also used by ThermalScalarModel"""
        default = {}
        default['frequency_choice']          = "end"
        default['frequency']                 = -1
        default['frequency_time']            = 1.
        default['format']                    = "ensight"
        default['time_dependency']           = 'fixed_mesh'
        default['options']                   = 'binary'
        default['repertory']                 = 'postprocessing'

        return default


    def __defaultWriterLabelAndId(self):
        """
        Private method.
        Return a default id and label for a new writer.
        """
        id_table = []
        for l in self.getWriterIdList():
            id_table.append(int(l))
        user_table = []
        for l in id_table:
            if l > 0:
                user_table.append(l)
        if user_table != []:
            next_id = max(user_table) +1
        else:
            next_id = 1
        next_label = 'writer('+ str(next_id)+')'
        n=next_id
        while next_label in self.getWriterLabelList():
            n=n+1
            next_label = 'writer('+ str(n)+')'
        return str(next_id), next_label


    def __updateWriterId(self):
        #"""
        #Private method.
        #Update suffixe number for writer label.
        #"""
        list = []
        n = 0
        for node in self.node_out.xmlGetNodeList('writer', 'label'):
            if int(node['id']) > 0 :
                n = n + 1
                if node['label'] == 'writer('+node['id']+')':
                    node['label'] ='writer('+str(n)+')'
                node['id'] = str(n)



    def addWriter(self):
        """Public method.
        Input a new user writer
        """

        i, l = self.__defaultWriterLabelAndId()
        if l not in self.getWriterIdList():
            self.node_out.xmlInitNode('writer', id = i,label = l)
        self.getWriterFrequencyChoice(i)
        self.getWriterFormat(i)
        self.getWriterRepertory(i)
        self.getWriterOptions(i)
        self.getWriterTimeDependency(i)

        return i



    def addDefaultWriter(self):
        """Public method.
        Input a new user writer
        """
        list_writer = self.getWriterIdList()
        if list_writer == []:
            node = self.node_out.xmlInitNode('writer', id = "-1", label = 'writer(-1)')
            node.xmlInitNode('frequency', choice = 'end')
            node.xmlInitNode('format', choice = 'ensight')
            node.xmlInitNode('repertory', choice = 'postprocessing')
            node.xmlInitNode('options', choice = 'binary')
            node.xmlInitNode('time_dependency', choice = 'fixed_mesh')



    def __deleteWriter(self, writer_id):
        """
        Private method.
        Delete writer.
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        node.xmlRemoveNode()
        self.__updateWriterId()


    def deleteWriter(self, writer_id):
        """
        Public method.
        Delete writer.
        """
        for w in writer_id:
            self.isInList(w, self.getWriterIdList())

        # First add the main scalar to delete
        list = writer_id
        # Delete all scalars
        for writer in list:
            self.__deleteWriter(writer)

        return list


    def getWriterIdList(self):
        """
        Return a list of writer id already defined
        """
        writer = []
        for node in self.node_out.xmlGetNodeList('writer', 'label'):
            writer.append(node['id'])
        return writer


    def getWriterLabelList(self):
        """
        Return a list of writer id already defined
        """
        writer = []
        for node in self.node_out.xmlGetNodeList('writer', 'label'):
            writer.append(node['label'])
        return writer


    def getWriterIdFromLabel(self, label):
        """
        Return the label of a writer
        """
        node = self.node_out.xmlGetNodeList('writer', 'label')
        for n in node:
            if n['label'] == label:
                writer_id = n['id']
        return writer_id


    def getWriterLabel(self, writer_id):
        """
        Return the label of a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        n = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        label = n['label']
        if label == None:
            label = __defaultWriterLabelAndId(id=writer_id)
            self.setWriterLabel(writer_id, label)
        return label


    def setWriterLabel(self, writer_id, label):
        """
        Set the label of a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        n = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n['label'] = label

    def getWriterFrequencyChoice(self, writer_id):
        """
        Return the choice of frequency output for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n= node.xmlGetNode('frequency_time')
        if n == None:
            n= node.xmlGetNode('frequency_formula')
            if n == None:
                n= node.xmlInitNode('frequency')
        frequency_choice = n['choice']
        if frequency_choice == None:
            frequency_choice = self.defaultWriterValues()['frequency_choice']
            self.setWriterFrequencyChoice(writer_id, frequency_choice)
        return frequency_choice


    def setWriterFrequencyChoice(self, writer_id, choice):
        """
        Set the choice of frequency output for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        self.isInList(choice, ('end', 'time', 'time_steps', 'second', 'formula'))
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        if choice == 'end' or choice == 'time' or choice == 'time_steps':
            n= node.xmlInitNode('frequency')
            n_bis = node.xmlGetNode('frequency_time')
            if n_bis != None:
                n_bis.xmlRemoveNode()
            n_ter = node.xmlGetNode('frequency_formula')
            if n_ter != None:
                n_ter.xmlRemoveNode()
        elif choice == 'second':
            n= node.xmlInitNode('frequency_time')
            n_bis = node.xmlGetNode('frequency')
            if n_bis != None:
                n_bis.xmlRemoveNode()
            n_ter = node.xmlGetNode('frequency_formula')
            if n_ter != None:
                n_ter.xmlRemoveNode()
        elif choice == 'formula':
            n= node.xmlInitNode('frequency_formula')
            n_bis = node.xmlGetNode('frequency_time')
            if n_bis != None:
                n_bis.xmlRemoveNode()
            n_ter = node.xmlGetNode('frequency')
            if n_ter != None:
                n_ter.xmlRemoveNode()
        n['choice'] = choice


    def getWriterFrequency(self, writer_id):
        """
        Return the frequency of a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        freq = node.xmlGetInt('frequency')
        if freq == None:
            freq = self.defaultWriterValues()['frequency']
            self.setWriterFrequency(writer_id, freq)
        return freq


    def setWriterFrequency(self, writer_id, freq):
        """
        Set the frequency of a writer
        """
        self.isInt(freq)
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        node.xmlSetData('frequency', freq)


    def getWriterFrequencyTime(self, writer_id):
        """
        Return the frequency of a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        freq = node.xmlGetDouble('frequency_time')
        if freq == None:
            freq = self.defaultWriterValues()['frequency_time']
            self.setWriterFrequencyTime(writer_id, freq)
        return freq


    def setWriterFrequencyTime(self, writer_id, freq):
        """
        Set the frequency of a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        node.xmlSetData('frequency_time', freq)


    def setWriterFrequencyFormula(self,writer_id, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        node.xmlSetData('frequency_formula', formula)


    def getWriterFrequencyFormula(self, writer_id):
        """
        Public method.
        Return the formula for a turbulent variable.
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        formula = node.xmlGetDouble('frequency')
        return formula


    def getWriterFormat(self, writer_id):
        """
        Return the format for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n= node.xmlInitNode('format')
        format = n['choice']
        if format == None:
            format = self.defaultWriterValues()['format']
            self.setWriterFormat(writer_id, format)
        return format


    def setWriterFormat(self, writer_id, format):
        """
        Set the format for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        self.isInList(format, ('ensight', 'med', 'cgns'))
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n= node.xmlInitNode('format')
        n['choice'] = format


    def getWriterRepertory(self, writer_id):
        """
        Return the repertory for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n= node.xmlInitNode('repertory')
        repertory = n['choice']
        if repertory == None:
            repertory = self.defaultWriterValues()['repertory']
            self.setWriterRepertory(writer_id, repertory)
        return repertory


    def setWriterRepertory(self, writer_id, repertory):
        """
        Set the repertory for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n= node.xmlInitNode('repertory')
        n['choice'] = repertory


    def getWriterOptions(self, writer_id):
        """
        Return the options for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n= node.xmlInitNode('options')
        options = n['choice']
        if options == None:
            options = self.defaultWriterValues()['options']
            self.setWriterOptions(writer_id, options)
        return options


    def setWriterOptions(self, writer_id, options):
        """
        Set the options for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n= node.xmlInitNode('options')
        n['choice'] = options


    def getWriterTimeDependency(self, writer_id):#-------> a réutiliser
        """
        Return the type of time dependency for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)

        n= node.xmlInitNode('time_dependency')
        choice = n['choice']
        if choice == None:
            choice = self.defaultWriterValues()['time_dependency']
            self.setWriterTimeDependency(writer_id, choice)
        return choice


    def setWriterTimeDependency(self, writer_id, choice):#-------> a réutiliser
        """
        Set the type of time dependency for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        self.isInList(choice, ('fixed_mesh', 'transient_coordinates', 'transient_connectivity'))
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n= node.xmlInitNode('time_dependency')
        n['choice'] = choice


    def defaultMeshValues(self):
        """Return the default values"""
        default = {}
        default['type']          = "cells"
        default['all_variables']          = "on"
        default['location']     = "all[]"
        default['time_dependendy']      = 'fixed_mesh'
        default['postprocessing_options'] = "binary"

        return default


    def getMeshIdList(self):
        """
        Return a list of mesh id already defined
        """
        mesh = []
        for node in self.node_out.xmlGetNodeList('mesh'):
            mesh.append(node["id"])
        return mesh


    def __defaultMeshLabelAndId(self):
        """
        Private method.
        Return a default id and label for a new mesh.
        """
        id_table = []
        for l in self.getMeshIdList():
            id_table.append(int(l))
        user_table = []
        for l in id_table:
            if l > 0:
                user_table.append(l)
        if user_table != []:
            next_id = max(user_table) +1
        else:
            next_id = 1
        next_label = 'mesh('+ str(next_id)+')'
        n=next_id
        while next_label in self.getMeshLabelList():
            n=n+1
            next_label = 'mesh('+ str(n)+')'
        return str(next_id), next_label


    def __updateMeshId(self):
        #"""
        #Private method.
        #Update suffixe number for mesh label.
        #"""
        list = []
        n = 0
        for node in self.node_out.xmlGetNodeList('mesh'):
            if int(node['id']) > 0 :
                n = n + 1
                if node['label'] == 'mesh('+node['id']+')':
                    node['label'] ='mesh('+str(n)+')'
                node['id'] = str(n)



    def addMesh(self):
        """Public method.
        Input a new user mesh
        """

        i, l = self.__defaultMeshLabelAndId()
        if l not in self.getMeshIdList():
            self.node_out.xmlInitNode('mesh', id = i,label = l)
        self.getMeshAllVariablesStatus(i)
        self.getMeshType(i)
        self.getMeshLocation(i)

        return i



    def addDefaultMesh(self):
        """Public method.
        Input default mesh
        """
        list_mesh = self.getMeshIdList()
        if list_mesh == []:
            node1 = self.node_out.xmlInitNode('mesh', id = "-1", label = 'fluid_domain', type = 'cells')
            node1.xmlInitNode('all_variables', status = 'on')
            node1.xmlInitNode('location')
            node1.xmlSetData('location','all[]')
            node1.xmlInitNode('writer', id = '-1')

            node2 = self.node_out.xmlInitNode('mesh', id = "-2", label = 'boundary_domain', type = 'boundary_faces')
            node2.xmlInitNode('all_variables', status = 'on')
            node2.xmlInitNode('location')
            node2.xmlSetData('location','all[]')
            node2.xmlInitNode('writer', id = '-1')


    def __deleteMesh(self, mesh_id):
        """
        Private method.
        Delete mesh.
        """
        self.isInList(mesh_id, self.getMeshIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        node.xmlRemoveNode()
        self.__updateMeshId()


    def deleteMesh(self, mesh_id):
        """
        Public method.
        Delete mesh.
        """
        for mesh in mesh_id:
            self.isInList(mesh, self.getMeshIdList())

        # First add the main scalar to delete
        list = mesh_id

        # Delete all scalars
        for mesh in list:
            self.__deleteMesh(mesh)

        return list


    def getMeshIdList(self):
        """
        Return a list of writer id already defined
        """
        mesh = []
        for node in self.node_out.xmlGetNodeList('mesh'):
            mesh.append(node['id'])
        return mesh



    def getMeshLabelList(self):
        """
        Return a list of mesh id already defined
        """
        mesh = []
        for node in self.node_out.xmlGetNodeList('mesh'):
            mesh.append(node['label'])
        return mesh


    def getMeshLabel(self, mesh_id):
        """
        Return the label of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        n = self.node_out.xmlGetNode('mesh', id = mesh_id)
        label = n['label']
        if label == None:
            label = __defaultMeshLabelAndId(id=mesh_id)
            self.setMeshLabel(mesh_id, label)
        return label


    def setMeshLabel(self, mesh_id, label):
        """
        Set the label of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        n = self.node_out.xmlGetNode('mesh', id = mesh_id)
        n['label'] = label


    def getMeshType(self, mesh_id):
        """
        Return the type of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        n = self.node_out.xmlGetNode('mesh', id = mesh_id)
        mesh_type = n['type']
        if mesh_type == None:
            mesh_type = self.defaultMeshValues()['type']
            self.setMeshType(mesh_id, mesh_type)
        return mesh_type


    def setMeshType(self, mesh_id, mesh_type):
        """
        Set the type of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        self.isInList(mesh_type, ('cells', 'interior_faces', 'boundary_faces'))
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        node['type'] = mesh_type


    def getMeshAllVariablesStatus(self, mesh_id):
        """
        Return the all_variables status of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        n= node.xmlInitNode('all_variables')
        status = n['status']
        if status == None:
            status = self.defaultMeshValues()['all_variables']
            self.setMeshAllVariablesStatus(mesh_id, status)
        return status


    def setMeshAllVariablesStatus(self, mesh_id, status):
        """
        Set the all_variables status of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        self.isOnOff(status)
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        n= node.xmlInitNode('all_variables')
        n['status'] = status

    def getMeshLocation(self, mesh_id):
        """
        Return the location of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        loc = node.xmlGetString('location')
        if loc == '':
            loc = self.defaultMeshValues()['location']
            self.setMeshLocation(mesh_id, loc)
        return loc


    def setMeshLocation(self, mesh_id, location):
        """
        Set the location of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        node.xmlSetData('location', location)


    def getAssociatedWriterIdList(self, mesh_id):
        """
        Return a list of associated writer to a mesh already defined
        """
        self.isInList(mesh_id, self.getMeshIdList())
        associated_writer = []
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        for n in node.xmlGetNodeList('writer'):
            associated_writer.append(n["id"])
        return associated_writer



    def addAssociatedWriter(self, mesh_id):
        """Public method.
        Input a new user associated writer to a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        n = self.node_out.xmlGetNode('mesh', id = mesh_id)
        writer_list = self.getWriterIdList()
        associated_writer_list = []
        for writer in writer_list:
            if writer not in self.getAssociatedWriterIdList(mesh_id):
                associated_writer_list.append(writer)
        writer_id = None
        if associated_writer_list:
            n.xmlInitNode('writer', id = associated_writer_list[0])
            writer_id = associated_writer_list[0]
        return writer_id


    def __deleteAssociatedWriter(self, mesh_id, writer_id):
        """
        Private method.
        Delete mesh.
        """
        self.isInList(mesh_id, self.getMeshIdList())
        self.isInList(writer_id, self.getAssociatedWriterIdList(mesh_id))
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        n = node.xmlGetNode('writer', id = writer_id)
        n.xmlRemoveNode()


    def deleteAssociatedWriter(self, mesh_id, writer_id):
        """
        Public method.
        Delete mesh.
        """
        self.isInList(mesh_id, self.getMeshIdList())
        self.isInList(writer_id, self.getAssociatedWriterIdList(mesh_id))

        self.__deleteAssociatedWriter(mesh_id, writer_id)

        return list


    def setAssociatedWriterChoice(self, mesh_id, writer_list):
        """
        Set the type of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        for w in writer_list :
            self.isInList(w, self.getWriterIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        for w in node.xmlGetNodeList('writer'):
            w.xmlRemoveNode()
        for w in writer_list:
            #print w, 'dans se associatedWriterchoice'
            node.xmlInitNode('writer', id = w)


    def getMonitoringPointType(self):
        """
        Return the type of output for printing listing
        """
        node = self.node_out.xmlGetNode('probe_recording_frequency')
        if node != None :
            return 'Frequency_h_x'
        val = self.getMonitoringPointFrequency()
        if val == -1 :
            return 'At the end'
        elif val == 1 :
            return 'At each step'
        else :
            return 'Frequency_h'


    def setMonitoringPointType(self, type):
        """
        Set the type of output for printing listing
        """
        self.isInList(type, ['None', 'At each step', 'Frequency_h', 'Frequency_h_x'])

        if type == 'Frequency_h_x' :
            childNode = self.node_out.xmlGetNode('probe_recording_frequency')
            if childNode != None :
                childNode.xmlRemoveNode()
        else :
            childNode = self.node_out.xmlGetNode('probe_recording_frequency_time')
            if childNode != None :
                childNode.xmlRemoveNode()


    def getMonitoringPointFrequency(self):
        """
        Return the frequency for recording probes
        """
        f = self.node_out.xmlGetInt('probe_recording_frequency')
        if f == None:
            f = self.defaultInitialValues()['probe_recording_frequency']
            self.setMonitoringPointFrequency(f)
        return f


    def setMonitoringPointFrequency(self, freq):
        """
        Set the frequency for recording probes
        """
        self.isInt(freq)
        self.node_out.xmlSetData('probe_recording_frequency', freq)


    def getMonitoringPointFrequencyTime(self):
        """
        Return the frequency for recording probes
        """
        f = self.node_out.xmlGetDouble('probe_recording_frequency_time')
        if f == None:
            f = self.defaultInitialValues()['probe_recording_frequency_time']
            self.setMonitoringPointFrequencyTime(f)
        return f


    def setMonitoringPointFrequencyTime(self, freq):
        """
        Set the frequency for recording probes
        """
        self.isFloat(freq)
        self.node_out.xmlSetData('probe_recording_frequency_time', freq)


    def getMonitoringPointFormat(self):
        """
        Return choice of format for post processing output file
        """
        node = self.node_out.xmlInitNode('probe_format', 'choice')
        choice = node['choice']
        if not choice:
            choice = self.defaultInitialValues()['probe_format']
            self.setMonitoringPointFormat(choice)
        return choice


    def setMonitoringPointFormat(self, choice):
        """
        Set choice of format for probes
        """
        self.isInList(choice, ('DAT', 'CSV'))
        node = self.node_out.xmlInitNode('probe_format', 'choice')
        node['choice'] = choice


    def addMonitoringPoint(self, x=0.0, y=0.0, z=0.0):
        """
        Public method.
        Add a new monitoring point.
        @type x: C{Float}
        @param x: first coordinate
        @type y: C{Float}
        @param y: second coordinate
        @type z: C{Float}
        @param z: third coordinate
        """
        self.isFloat(x)
        self.isFloat(y)
        self.isFloat(z)
        num = str(self.getNumberOfMonitoringPoints() + 1)
        status="on"
        node = self.node_out.xmlInitNode('probe', name=num, status=status)
        for coord, val in [('probe_x', x), ('probe_y', y), ('probe_z', z)]:
            self.__setCoordinates(num, coord, val)


    def replaceMonitoringPointCoordinates(self, name, x=0.0, y=0.0, z=0.0):
        """
        Public method.
        Change the coordinates of a monitoring point.
        @type name: C{String}
        @param name: identifier of the monitoring point
        @type x: C{Float}
        @param x: first new coordinate
        @type y: C{Float}
        @param y: second new coordinate
        @type z: C{Float}
        @param z: third new coordinate
        """
        self.isFloat(x)
        self.isFloat(y)
        self.isFloat(z)
        self.isStr(name)
        self.isGreater(float(name), 0.0)
        self.isLowerOrEqual(float(name), self.getNumberOfMonitoringPoints())

        for coord, val in [('probe_x', x), ('probe_y', y), ('probe_z', z)]:
            self.__setCoordinates(name, coord, val)


    def deleteMonitoringPoints(self, list):
        """
        Public method.
        Conveniant method for the view. Delete a list of monitoring points.
        @type list: C{List} of C{Int}
        @param list: list of identifier of monitoring points to delete
        """
        list.sort()
        r = len(list)
        for n in range(r):
            name = str(list[n])
            self.deleteMonitoringPoint(name)
            for i in range(n, r):
                list[i] = list[i] - 1


    def deleteMonitoringPoint(self, num):
        """
        Public method.
        Delete a single monitoring point.
        @type num: C{String}
        @param num: identifier of the monitoring point
        """
        self.isStr(num)
        self.isGreater(float(num), 0.0)
        self.isLowerOrEqual(float(num), self.getNumberOfMonitoringPoints())

        # delete the node of the monitoring point

        node = self.node_out.xmlGetNode('probe', name=num)
        if node:
            node.xmlRemoveNode()
            self.case.xmlRemoveChild('probe_recording', name=num)

            # renumerotation of all monitoring points

            for p in range(int(num)+1, self.getNumberOfMonitoringPoints()+2):
                probe = self.node_out.xmlGetNode('probe', name=p)
                probe['name'] = p - 1
                for probe_recording in self.case.xmlGetNodeList('probe_recording', name=p):
                    probe_recording['name'] = p - 1

            # update the attribute "choice" of the probes markup for variables

            from Pages.OutputVolumicVariablesModel import OutputVolumicVariablesModel
            listNodeVolum = OutputVolumicVariablesModel(self.case).listNodeVolum
            del OutputVolumicVariablesModel
            for nodeList in listNodeVolum:
                for node in nodeList:
                    n = node.xmlGetChildNode('probes')
                    if n:
                        nlist = n.xmlGetChildNodeList('probe_recording')
                        if not nlist:
                            n.xmlRemoveNode()
                        else:
                            n['choice']= str(len(nlist))


    def getMonitoringPointCoordinates(self, name):
        """
        Public method.
        @type name: C{String}
        @param name: identifier of the monitoring point
        @return: coordinates X, Y, Z for the monitoring point I{name}
        @rtype: C{List} of C{Float}
        """
        self.isStr(name)
        self.isGreater(float(name), 0.0)
        self.isLowerOrEqual(float(name), self.getNumberOfMonitoringPoints())
        X = self.__getCoordinates(name, 'probe_x')
        Y = self.__getCoordinates(name, 'probe_y')
        Z = self.__getCoordinates(name, 'probe_z')
        return X, Y, Z


    def getNumberOfMonitoringPoints(self):
        """
        Public method.
        @return: number of monitoring points already defined.
        @rtype: C{Int}
        """
        return len(self.node_out.xmlGetNodeList('probe'))


#-------------------------------------------------------------------------------
# OutputControlModel test Class
#-------------------------------------------------------------------------------

class OutputControlModelTestCase(ModelTest):
    """
    """
    def checkOutputControlModeInstantiation(self):
        """Check whether the OutputControlModel class could be instantiated"""
        model = None
        model = OutputControlModel(self.case)
        assert model != None, 'Could not instantiate OutputControlModel'

    def checkSetandGetListingFrequency(self):
        """Check whether the frequency of output listing could be set and get"""
        model = OutputControlModel(self.case)
        model.setListingFrequency(12)
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <listing_printing_frequency>12</listing_printing_frequency>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
                    'Could not set frequency for listing output control model'
        assert model.getListingFrequency() == 12, \
                    'Could not get frequency for listing output control model'

    def checkSetandGetPostprocessingFrequency(self):
        """Check whether the frequency of post processing could be set and get"""
        model = OutputControlModel(self.case)
        model.setPostprocessingFrequency(13)
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <postprocessing_frequency>13</postprocessing_frequency>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
                    'Could not set frequency for post processing output control model'
        assert model.getPostprocessingFrequency() == 13, \
                    'Could not get frequency for post processing output control model'

    def checkSetandGetFluidDomainPostProStatus(self):
        """Check whether the status of post processing for fluid domain could be set and get"""
        model = OutputControlModel(self.case)
        model.setFluidDomainPostProStatus('off')
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <fluid_domain status="off"/>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
            'Could not set status of post processing for fluid domain for output control model'
        assert model.getFluidDomainPostProStatus() == 'off', \
            'Could not get status of post processing for fluid domain for output control model'

    def checkSetandGetDomainBoundaryPostProStatus(self):
        """
        Check whether the status of post processing for domain
        boundary could be set and get
        """
        model = OutputControlModel(self.case)
        model.setDomainBoundaryPostProStatus('on')
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <domain_boundary status="on"/>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
        'Could not set status of post processing for domain boundary for output control model'
        assert model.getDomainBoundaryPostProStatus() == 'on', \
        'Could not get status of post processing for domain boundary for output control model'

    def checkSetandGetTypePostMeshes(self):
        """Check whether the type of mesh's post processing could be set and get"""
        model = OutputControlModel(self.case)
        model.setTypePostMeshes('2')
        doc = '''<output>
                    <postprocessing_mesh_options choice="2"/>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
        'Could not set type of post processing for mesh in output control model'
        assert model.getTypePostMeshes() == '2', \
        'Could not get type of post processing for mesh in output control model'

    def checkSetandGetPostProFormat(self):
        """Check whether the format for post processing could be set and get"""
        model = OutputControlModel(self.case)
        model.setPostProFormat('MED')
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <postprocessing_format choice="MED"/>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
        'Could not set format of post processing in output control model'
        assert model.getPostProFormat() == 'MED', \
        'Could not get format of post processing in output control model'

    def checkSetandGetPostProOptionsFormat(self):
        """
        Check whether the options of format for post processing could be set and get
        """
        model = OutputControlModel(self.case)
        model.setPostProOptionsFormat('big_endian,divide_polyhedra')
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <postprocessing_options choice="big_endian,divide_polyhedra"/>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
        'Could not set the options of format for post processing in output control model'
        assert model.getPostProOptionsFormat() == 'big_endian,divide_polyhedra', \
        'Could not get the options of format for post processing in output control model'

    def checkSetandGetMonitoringPointFrequency(self):
        """
        Check whether the frequency of monitoring point could be set and get
        """
        model = OutputControlModel(self.case)
        model.setMonitoringPointFrequency(15)
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <probe_recording_frequency>15</probe_recording_frequency>
                  </output>'''
        assert model.node_out== self.xmlNodeFromString(doc),\
        'Could not set the frequency of monitoring point in output control model'
        assert model.getMonitoringPointFrequency() == 15,\
        'Could not get the frequency of monitoring point in output control model'

    def checkAddMonitoringPoint(self):
        """
        Check whether monitoring point could be added
        """
        model = OutputControlModel(self.case)
        model.addMonitoringPoint(11.1, 22.2, 33.3)
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <probe name="1" status="on">
                        <probe_x>11.1</probe_x>
                        <probe_y>22.2</probe_y>
                        <probe_z>33.3</probe_z>
                    </probe>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc),\
        'Could not add monitoring point in output control model'
        assert model.getMonitoringPointCoordinates("1") == (11.1, 22.2, 33.3),\
        'Could not get monitoring point in output control model'

    def checkReplaceMonitoringPoint(self):
        """
        Check whether monitoring point could be replaced
        """
        model = OutputControlModel(self.case)
        model.addMonitoringPoint(11.1, 22.2, 33.3)
        model.addMonitoringPoint(5, 5.1, 5.21)
        model.replaceMonitoringPointCoordinates("2",5., 6, 7.)
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <probe name="1" status="on">
                        <probe_x>11.1</probe_x>
                        <probe_y>22.2</probe_y>
                        <probe_z>33.3</probe_z>
                    </probe>
                    <probe name="2" status="on">
                        <probe_x>5</probe_x>
                        <probe_y>6</probe_y>
                        <probe_z>7</probe_z>
                    </probe>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc),\
        'Could not replace monitoring point in output control model'


    def checkDeleteMonitoringPoint(self):
        """
        Check whether monitoring point could be deleted
        """
        model = OutputControlModel(self.case)
        model.addMonitoringPoint(11.1, 22.2, 33.3)
        model.addMonitoringPoint(5, 5.1, 5.21)
        model.addMonitoringPoint(9.,8.,7.)
        model.deleteMonitoringPoint("2")
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <probe name="1" status="on">
                        <probe_x>11.1</probe_x>
                        <probe_y>22.2</probe_y>
                        <probe_z>33.3</probe_z>
                    </probe>
                    <probe name="2" status="on">
                        <probe_x>9</probe_x>
                        <probe_y>8</probe_y>
                        <probe_z>7</probe_z>
                    </probe>
                </output>'''
        assert model.node_out== self.xmlNodeFromString(doc),\
        'Could not delete monitoring point in output control model'


def suite():
    testSuite = unittest.makeSuite(OutputControlModelTestCase, "check")
    return testSuite

def runTest():
    print("OutputControlModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
