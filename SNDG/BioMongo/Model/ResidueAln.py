'''
Created on Mar 22, 2016

@author: eze
'''
from mongoengine.fields import IntField

from SNDG.BioMongo.Model.Alignment import SimpleAlignment, AlnPos


class AlnPosStructural(AlnPos):
    '''
    '''

    def __init__(self, pos, qpos, qaa, hpos, haa, q_resid, h_resid, q_offset=0, h_offset=0):
        AlnPos.__init__(self, pos, qpos, qaa, hpos, haa, q_offset=0, h_offset=0)
        self.q_resid = q_resid
        self.h_resid = h_resid

    def __str__(self):
        return "AlnPos(pos=%i,q_resid=%i,h_resid=%i)" % (self.pos, self.q_resid, self.h_resid)


class ResidueAln(SimpleAlignment):
    '''
    classdocs
    '''
    query_res_start = IntField(default=1)
    query_res_end = IntField()
    hit_res_start = IntField(default=1)
    hit_res_end = IntField()

    def _init_pos_map(self):
        SimpleAlignment._init_pos_map(self)
        old_pos = self.positions
        self.positions = []

        for j in range(self.aln_len()):
            aln_pos = AlnPosStructural(j,
                                       old_pos[j]._q_pos,
                                       old_pos[j].qaa,
                                       old_pos[j]._h_pos,
                                       old_pos[j].haa,
                                       self.query_res_start + j - old_pos[j].q_offset if not self.aln_query.is_gap(
                                           j) else None,
                                       # + old_pos[j]._q_pos,
                                       self.hit_res_start + j - old_pos[j].h_offset if not self.aln_hit.is_gap(
                                           j) else None,
                                       # old_pos[j]._h_pos  if not self.aln_hit.is_gap(j) else None,
                                       old_pos[j].q_offset,
                                       old_pos[j].h_offset)
            self.positions.append(aln_pos)

    def aln_pos_from_query_pos(self, query_pos):
        return [x for x in self.positions if x._q_pos == query_pos][0]

    def aln_pos_from_q_resid(self, q_resid):
        return [x for x in self.positions if x.q_resid == q_resid][0]

    def aln_pos_from_h_resid(self, h_resid):
        return [x for x in self.positions if x.h_resid == h_resid][0]
