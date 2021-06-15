# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 23:46:33 2021

@author: Admin
"""

import numpy as np

class Vershamov_Tenengoltz:
    def __init__(self, k: int, q: int, a = 0):
            assert q >= 2   #checking valid input
            assert k > 0
            
            self.q = q
            self.k = k
            self.n = self._compute_n_for_k()
            self.a = a
            
            if q==2:
                self.m = self.n + 1
                self._generate_non_dyadic_positions()
            else:
                self.m = self.n
                self.t = np.ceil(np.log2(self.n)).astype(np.int64)
                assert 0 <= self.a < self.m
                assert 0 <= self.q
                self._generate_maps()
                self._generate_q_ary_step_1_positions_pos()
            
    
    def encode(self, message):
        message = np.array(message, dtype = np.int64)
        assert message.ndim == 1
        assert message.size == self.k
        if self.q==2:
            codeword=self.encode_binary(message)
        else:
            codeword=self._encode_q_ary(message)
        return codeword
    
    def decode(self, codeword):
        codeword = np.array(codeword, dtype=np.int64)
        assert codeword.ndim == 1
        codeword_len = codeword.size
        if (codeword_len <= self.n + 1 and codeword_len >= self.n - 1):
            if(np.max(codeword) < self.q and np.min(codeword) > -1):
                if(self.q == 2):
                    if codeword_len != self.n:
                        codeword = self._fix_binary(codeword)
                    if not self._is_codeword_with_syn_a(codeword):
                        return None
                    return self._decode_binary(codeword)
                else:
                    if codeword_len != self.n:
                        codeword = self._fix_q_array(codeword)
                    if not self._is_codeword_with_syn_a(codeword):
                        return None
                    return self._decode_q_ary(codeword)
        return None     #more then 1 substitution or insertion, or unvalid values in the codeword
    
    def encode_binary(self,message):
        codeword = np.zeros(self.n, dtype = np.int64)
        codeword[self.non_dyadic_positions-1] = message
        deficiency = self._compute_deficiency_binary(codeword)
        if deficiency != 0:
            for pos in range(len(self.dyadic_positions)):
                codeword[self.dyadic_positions[pos] - 1] = deficiency%2
                deficiency = np.floor(deficiency/2)
                if deficiency == 0:
                    break
        assert self._is_codeword_with_syn_a(codeword)
        return codeword
    
    
    def _encode_q_ary(self, message):

        codeword = np.zeros(self.n, dtype = np.int64)
        # step 1
        read_bits = np.floor(np.max([self.n-3*self.t+3,0])*np.log2(self.q)).astype(np.int64)
        if read_bits > 0:
            codeword[self.q_ary_step_1_positions] = self._convert_base(message[:read_bits],
                    2, self.q, out_size = self.q_ary_step_1_positions.size)
        
        
        # step 2
        bits_per_tuple_step_2 = np.floor(2*np.log2(self.q-1)).astype(np.int64)
        for i in range(3, self.t):
            map_for_pairs_index = self.convert_from_q_array(message[read_bits:read_bits+bits_per_tuple_step_2], 2)
            if 2**i == self.n-1:
                num_bits_special_case = np.floor(np.log2(self.q-1)).astype(np.int64)
                codeword[2**i-1] = self.convert_from_q_array(message[read_bits:read_bits+num_bits_special_case], 2)+1
                read_bits += num_bits_special_case
                break
            codeword[2**i-1] = self.map_for_pairs_r[map_for_pairs_index]
            codeword[2**i+1] = self.map_for_pairs_l[map_for_pairs_index]
            read_bits += bits_per_tuple_step_2
        
        ##special case (s_3, s_5)
        codeword[3] = self.q - 1
        if self.q == 3:
            codeword[5] = 2
        else:
            bits_for_s_5 = np.floor(np.log2(self.q-1)).astype(np.int64)
            map_for_s_5_index = self.convert_from_q_array(message[read_bits:read_bits+bits_for_s_5],2)
            codeword[5] = self.map_for_s_5[map_for_s_5_index]
            read_bits += bits_for_s_5
        assert read_bits == self.k
        
        
        # step 3
        auxiliary_sequence = self._convert_codeword_to_auxiliary_sequence(codeword)
        for i in range(2, self.t):
            if 2**i == self.n-1:
                break
            auxiliary_sequence[2**i+1-1] = (codeword[2**i+1] >= codeword[2**i-1])
        auxiliary_sequence[3-1] = 1
        
        
        
        # step 4
        for i in range(self.t):
            auxiliary_sequence[2**i-1] = 0
        deficiency = self._compute_deficiency_binary(auxiliary_sequence)
        if deficiency != 0:
            for i in range(self.t):
                auxiliary_sequence[2**i - 1] = deficiency%2
                deficiency = np.floor(deficiency/2)
                if deficiency == 0:
                    break
        
        # step 5
        for i in range(2,self.t):
            pos = 2**i
            if auxiliary_sequence[pos-1] == 0:
                codeword[pos] = codeword[pos-1]-1
            else:
                codeword[pos] = codeword[pos-1]
        
        
        # step 6
        w = np.mod(0-np.sum(codeword[3:]),self.q)
        if self.q == 3:
            if auxiliary_sequence[1-1] == 1 and auxiliary_sequence[2-1] == 1:
                codeword[2], codeword[1], codeword[0] = 2, 2, np.mod(w-4, 3)
            elif auxiliary_sequence[1-1] == 1 and auxiliary_sequence[2-1] == 0:
                codeword[2], codeword[1], codeword[0] = 1, 2, w
            elif auxiliary_sequence[1-1] == 0 and auxiliary_sequence[2-1] == 1:
                codeword[2] = 2
                if w == 1:
                    codeword[1], codeword[0] = 0, 2
                elif w == 0:
                    codeword[1], codeword[0] = 0, 1
                else:
                    codeword[1], codeword[0] = 1, 2
            else:
                auxiliary_sequence[1-1], auxiliary_sequence[2-1], auxiliary_sequence[3-1] = 1, 1, 0
                codeword[3] = 1
                codeword[4] = 0 if (auxiliary_sequence[4-1] == 0) else 1
                codeword[2], codeword[1] = 2, 2
                codeword[0] = np.mod(0-np.sum(codeword[1:]), 3)
        else:
            if w == 1:
                x, y, z = 0, 2, self.q - 1
            elif w == 2:
                x, y, z = 1, 2, self.q - 1
            else:
                x, y, z = 0, 1, np.mod(w-1, self.q)
            if auxiliary_sequence[1-1] == 0 and auxiliary_sequence[2-1] == 0:
                codeword[0], codeword[1], codeword[2] = z, y, x
            elif  auxiliary_sequence[1-1] == 0 and auxiliary_sequence[2-1] == 1:
                codeword[0], codeword[1], codeword[2] = z, x, y
            elif  auxiliary_sequence[1-1] == 1 and auxiliary_sequence[2-1] == 0:
                codeword[0], codeword[1], codeword[2] = x, z, y
            else:
                codeword[0], codeword[1], codeword[2] = x, y, z
        assert np.array_equal(auxiliary_sequence, self._convert_codeword_to_auxiliary_sequence(codeword))
        assert self._is_codeword_with_syn_a(codeword)
        return codeword



    def _decode_binary(self, codeword):

        return codeword[self.non_dyadic_positions-1]

    def _decode_q_ary(self, codeword):

        message = np.zeros(self.k, dtype = np.int64)
        # step 1
        read_bits = np.floor(np.max([self.n-3*self.t+3,0])*np.log2(self.q)).astype(np.int64)
        if read_bits > 0:
            step_1_bits = self._convert_base(codeword[self.q_ary_step_1_positions], self.q, 2, read_bits)
            if step_1_bits is None: 
                return None
            message[:read_bits] = step_1_bits

        # step 2
        bits_per_tuple_step_2 = np.floor(2*np.log2(self.q-1)).astype(np.int64)
        for i in range(3, self.t):
            if 2**i == self.n-1:
                num_bits_special_case = np.floor(np.log2(self.q-1)).astype(np.int64)
                if codeword[2**i-1] == 0:
                    return None
                message[read_bits:read_bits+num_bits_special_case] = self.convert_to_q_array(codeword[2**i-1]-1, 2, num_bits_special_case)
                read_bits += num_bits_special_case
                break
            if (codeword[2**i-1], codeword[2**i+1]) in self.map_for_pairs_rev:
                message[read_bits:read_bits+bits_per_tuple_step_2] = self.convert_to_q_array(self.map_for_pairs_rev[(codeword[2**i-1], codeword[2**i+1])], 2, bits_per_tuple_step_2)
            else:
                return None
            read_bits += bits_per_tuple_step_2

            ##special case q = 3
        if self.q == 3:
            if codeword[5] != 2:
                return None
        else:
            if codeword[3] != self.q - 1:
                return None
            bits_for_s_5 = np.floor(np.log2(self.q-1)).astype(np.int64)
            if codeword[5] not in self.map_for_s_5_rev:
                return None
            else:
                message[read_bits:read_bits+bits_for_s_5] = self.convert_to_q_array(self.map_for_s_5_rev[codeword[5]], 2, bits_for_s_5)
            read_bits += bits_for_s_5
        assert read_bits == self.k
        return message
    
    def _fix_q_array(self, codeword):
        auxiliary_sequence = self._convert_codeword_to_auxiliary_sequence(codeword)
        auxiliary_sequence_corrected = self._fix_binary(auxiliary_sequence)
        if auxiliary_sequence_corrected is None or self._compute_deficiency_binary(auxiliary_sequence_corrected) != 0:
            return None
        codeword_decoded = np.zeros(self.n, dtype=np.int64)
        if auxiliary_sequence.size == self.n-2:
            # deletion
            error_symbol = np.mod(0-np.sum(codeword),self.q)
            if np.array_equal(auxiliary_sequence, auxiliary_sequence_corrected[:-1]):
                diff_pos = self.n-2
            else:
                for diff_pos in range(self.n-2):
                    if auxiliary_sequence[diff_pos] != auxiliary_sequence_corrected[diff_pos]:
                        break

            del_pos_found = False
            for del_pos in reversed(range(diff_pos + 2)):
                if del_pos == 0:
                    if auxiliary_sequence_corrected[0] == (codeword[0] >= error_symbol):
                        del_pos_found = True
                        break
                elif del_pos == self.n-1:
                    if (auxiliary_sequence_corrected[self.n-2] == (error_symbol >= codeword[self.n-2])):
                        del_pos_found = True
                        break
                else:
                    if (auxiliary_sequence_corrected[del_pos-1] == (error_symbol >= codeword[del_pos-1])) \
                    and (auxiliary_sequence_corrected[del_pos+1-1] == (codeword[del_pos] >= error_symbol)):
                        del_pos_found = True
                        break
            if del_pos_found:
                codeword_decoded[:del_pos] = codeword[:del_pos]
                codeword_decoded[del_pos] = error_symbol
                codeword_decoded[del_pos+1:] = codeword[del_pos:]
            else:
                codeword_decoded = None
        else:
            #insertion
            error_symbol = np.mod(np.sum(codeword),self.q)
            if np.array_equal(auxiliary_sequence[:-1], auxiliary_sequence_corrected):
                diff_pos = self.n-1
            else:
                for diff_pos in range(self.n):
                    if auxiliary_sequence[diff_pos] != auxiliary_sequence_corrected[diff_pos]:
                        break
            ins_pos_found = False
            for ins_pos in reversed(range(diff_pos + 2)):
                if ins_pos == 0 or ins_pos == self.n:
                    if (codeword[ins_pos] == error_symbol):
                        ins_pos_found = True
                        break
                else:
                    if (codeword[ins_pos] == error_symbol) and \
                        (auxiliary_sequence_corrected[ins_pos-1] == (codeword[ins_pos+1] >= codeword[ins_pos-1])):
                        ins_pos_found = True
                        break
            if ins_pos_found:
                codeword_decoded[:ins_pos] = codeword[:ins_pos]
                codeword_decoded[ins_pos:] = codeword[ins_pos+1:]
            else:
                codeword_decoded = None
        if codeword_decoded is not None and self._compute_deficiency_q_ary(codeword_decoded) == (0,0):
            return codeword_decoded
        else:
            return None
        
        
    def _fix_binary(self, codeword):
    
        s = self._compute_deficiency_binary(codeword)
        w = np.sum(codeword)
        codeword_decoded = np.zeros(self.n, dtype=np.int64)
        if codeword.size == self.n-1:
            # deletion
            if s == 0:
                codeword_decoded[:-1] = codeword
            elif s <= w:
                num_ones_seen = 0
                for i in reversed(range(self.n-1)):
                    if codeword[i] == 1:
                        num_ones_seen += 1
                        if num_ones_seen == s:
                            codeword_decoded[:i] = codeword[:i]
                            codeword_decoded[i+1:] = codeword[i:]
                            break
            else:
                num_zeros_seen = 0
                if s-w-1 == 0:
                    codeword_decoded[0] = 1
                    codeword_decoded[1:] = codeword
                else:
                    success = False
                    for i in range(self.n-1):
                        if codeword[i] == 0:
                            num_zeros_seen += 1
                            if num_zeros_seen == s-w-1:
                                codeword_decoded[:i+1] = codeword[:i+1]
                                codeword_decoded[i+1] = 1
                                codeword_decoded[i+2:] = codeword[i+1:]
                                success = True
                                break
                    if not success:
                        codeword_decoded = None
        else:
            # insertion
            if s == self.m-self.n-1 or s == 0:
                codeword_decoded = codeword[:-1]
            elif s == self.m-w:
                codeword_decoded = codeword[1:]
            elif s > self.m-w:
                num_ones_seen = 0
                success = False
                for i in reversed(range(2,self.n+1)):
                    if codeword[i] == 1:
                        num_ones_seen += 1
                        if num_ones_seen == self.m-s:
                            if codeword[i-1] == 0:
                                codeword_decoded[:i-1] = codeword[:i-1]
                                codeword_decoded[i-1:] = codeword[i:]
                                success = True
                            else:
                                pass
                            break
                if not success:
                    codeword_decoded = None
            else:
                num_zeros_seen = 0
                success = False
                for i in range(self.n-1):
                    if codeword[i] == 0:
                        num_zeros_seen += 1
                        if num_zeros_seen == self.m-w-s:
                            if codeword[i+1] == 1:
                                codeword_decoded[:i+1] = codeword[:i+1]
                                codeword_decoded[i+1:] = codeword[i+2:]
                                success = True
                            else:
                                pass
                            break
                if not success:
                    codeword_decoded = None
        return codeword_decoded

    def _compute_deficiency_q_ary(self, codeword):
        auxiliary_sequence = self._convert_codeword_to_auxiliary_sequence(codeword)
        return (self._compute_deficiency_binary(auxiliary_sequence), np.mod(0-np.sum(codeword),self.q))#always decode to a VT group of sum = 0
   
    
    def _is_codeword_with_syn_a(self, codeword):
        if codeword is None or codeword.size != self.n:
            return False
        if self.q == 2:
            return (self._compute_deficiency_binary(codeword) == 0)
        else:
            return (self._compute_deficiency_q_ary(codeword) == (0,0))
    
    def _generate_non_dyadic_positions(self):
        t = np.ceil(np.log2(self.n+1)).astype(np.int64)
        self.dyadic_positions = np.zeros(self.n-self.k, dtype=np.int64)
        for i in range(t):
            self.dyadic_positions[i] = 2**i
        self.non_dyadic_positions =  np.setdiff1d(np.arange(1,self.n+1), self.dyadic_positions)
        return

    def _compute_n_for_k(self):
        n = self.k + np.ceil(np.log2(self.k+1)).astype(np.int64)
        while True:
            if self._find_k(n) >= self.k:
                break
            n += 1
        return n

    def _find_k(self,n):
        return n - np.ceil(np.log2(n+1)).astype(np.int64)

    
    def convert_from_q_array(self,q_ary_array, q):
        num = 0
        for i in q_ary_array:
            i = i.item()
            num = q*num + i
        return num

    def convert_to_q_array(self,num, q, array_size = None):

        q_base_array = []
        while num > 0:
            q_base_array.append(num%q)
            num = np.floor(num/q)
        q_base_array.reverse()
        q_base_array = np.array(q_base_array, dtype = np.int64)
        if array_size == None:
            return q_base_array
        if q_base_array.size > array_size:
            return None
        else:
            return np.pad(q_base_array, (array_size - q_base_array.size,0))
    
    def _convert_base(self,in_array, in_base, out_base, out_size = None):
        # convert array in in_base to num and then num to array in out_base
        num = self.convert_from_q_array(in_array, in_base)
        return self.convert_to_q_array(num, out_base, out_size)

    def _compute_deficiency_binary(self, codeword):
        codeword_len = codeword.size
        return np.mod(self.a - np.sum((1+np.arange(codeword_len))*codeword),self.m)
    
    def create_deletion(self,codeword):
        return np.delete(codeword, np.random.randint(0,len(codeword)))
    
    def _convert_codeword_to_auxiliary_sequence(self,codeword):
        return (codeword[1:] >= codeword[:-1]).astype(np.int64)
    
    def _generate_maps(self):
      
      map_for_pairs_size = 2**(np.floor(2*np.log2(self.q-1)).astype(np.int64))
      self.map_for_pairs_l = np.zeros(map_for_pairs_size, dtype = np.int64)
      self.map_for_pairs_r = np.zeros(map_for_pairs_size, dtype = np.int64)
      counter = 0
      for r in range(self.q):
          if counter == map_for_pairs_size:
              break
          if r == 0:        # r != 0 and 
              continue
          for l in range(self.q):
              if counter == map_for_pairs_size:
                  break
              if l == r-1:      # l != r-1
                  continue
              self.map_for_pairs_l[counter] = l
              self.map_for_pairs_r[counter] = r
              counter += 1

      self.map_for_pairs_rev = {(self.map_for_pairs_r[i],self.map_for_pairs_l[i]): i for i in range(map_for_pairs_size)}

      if self.q != 3: # 3 is special base
          map_for_s_5_size = 2**(np.floor(np.log2(self.q-1)).astype(np.int64))
          self.map_for_s_5 = np.zeros(map_for_s_5_size, dtype = np.int64)
          for i in range(map_for_s_5_size):
              self.map_for_s_5[i] = i
              if i == self.q-2:
                  self.map_for_s_5[i] = self.q-1
          
          self.map_for_s_5_rev = {self.map_for_s_5[i]: i for i in range(map_for_s_5_size)}
      return
    
    def _generate_q_ary_step_1_positions_pos(self):
      non_q_ary_step_1_positions_pos = [0, 1, 2]
      for i in range(2,self.t):
          non_q_ary_step_1_positions_pos += [(2**i)-1, 2**i, (2**i)+1]
      non_q_ary_step_1_positions_pos = np.array(non_q_ary_step_1_positions_pos, dtype = np.int64)
      self.q_ary_step_1_positions =  np.setdiff1d(np.arange(self.n), non_q_ary_step_1_positions_pos)

