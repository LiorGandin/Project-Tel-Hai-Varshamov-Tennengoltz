# -*- coding: utf-8 -*-
"""
Created on Mon May 17 12:48:32 2021

@author: Admin
"""

import numpy as np

class PreCodingDataBreakDown:
    
    def split(self,word):
        return [char for char in word]
    
    def PreCodingDataBreakDownSpliter(self,a_string,chunk_size,q):
        a_byte_array = bytearray(a_string, "utf8")
        byte_list = []
        for byte in a_byte_array:
            q_representation = np.base_repr(byte,q)
            byte_list.append(q_representation)
            
        splited_array=[]
        for i in range(0,len(byte_list)):
            for j in range(0,len(byte_list[i]),chunk_size):
                end_of_chunk=j+chunk_size
                if j+chunk_size>len(byte_list[i]):
                    end_of_chunk=len(byte_list[i])
                    tmp=list(np.zeros(chunk_size-(end_of_chunk-j)))
                    tmp=tmp+self.split(byte_list[i][j:end_of_chunk])
                else:
                    tmp=self.split(byte_list[i][j:end_of_chunk])
                tmp = [int(i) for i in tmp]
                splited_array=splited_array+tmp
        return splited_array