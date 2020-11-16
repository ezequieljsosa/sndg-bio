'''
Created on Aug 14, 2015

@author: eze
'''
from mongoengine.document import EmbeddedDocument
from mongoengine.fields import IntField, FloatField, EmbeddedDocumentField, \
    StringField

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def map_alignment(query_aln, subject_aln, gap_code="-"):
    '''
    '''
    pos_map = {}
    i = 0
    query_gaps = 0
    subject_gaps = 0
    for j in range(len(query_aln)):
        #        assert query_aln[j] != "-", "No puede ser que tenga un agujero el RV"
        if query_aln[j] == gap_code:
            query_gaps = query_gaps + 1
        query_pos = j - query_gaps

        if subject_aln[j] == gap_code:
            pos_map[query_pos] = None
        else:
            pos_map[query_pos] = i
            i = i + 1
    return pos_map


class AlnPos(object):
    '''
    '''

    def __init__(self, pos, qpos, qaa, hpos, haa, q_offset=0, h_offset=0):
        '''
        Constructor
        '''
        self.pos = pos
        self._q_pos = qpos
        self.qaa = qaa
        self._h_pos = hpos
        self.haa = haa
        self.q_offset = q_offset
        self.h_offset = h_offset

    def q_pos(self):
        return self._q_pos + self.q_offset

    def h_pos(self):
        return self._h_pos + self.h_offset

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "AlnPos(" + str(self.pos) + " " + str(self._q_pos) + ":" + str(self.qaa) + " " + str(
            self._h_pos) + ":" + str(self.haa) + ")"


class Alignment(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.aln_lines = []


class AlnLine(EmbeddedDocument):
    meta = {'allow_inheritance': True}

    name = StringField(required=True)
    seq = StringField(required=True)
    start = IntField(default=0)
    end = IntField()
    strand = IntField(default=1)

    def gaps(self):
        return self.seq.count("-")

    def original_seq(self):
        return self.seq.replace("-", "")

    def is_gap(self, pos):
        return self.seq[pos] == "-"

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "AlnLine(" + self.name + ":" + str(self.start) + "-" + str(self.end) + ")"


class SimpleAlignment(EmbeddedDocument):
    '''
    classdocs
    '''
    meta = {'allow_inheritance': True}

    evalue = FloatField()
    identity = FloatField()
    query_coverage = FloatField()
    hit_coverage = FloatField()
    aln_query = EmbeddedDocumentField(AlnLine, required=True)
    aln_hit = EmbeddedDocumentField(AlnLine, required=True)
    aln_mid = StringField()
    aln_cd = StringField()
    aln_pp = StringField()

    def __init__(self, **kwargs):
        '''
        Constructor
        '''
        super(EmbeddedDocument, self).__init__(**kwargs)
        self._instance = None
        self._init_pos_map()
        if not self.identity:
            self.update_identity()

    def _init_pos_map(self):
        self.positions = []

        q_start = self.aln_query.start if self.aln_query.strand == 1 else self.aln_query.end
        h_start = self.aln_hit.start if self.aln_hit.strand == 1 else self.aln_hit.end

        query_gaps = 0
        hit_gaps = 0
        for j in range(self.aln_len()):
            if self.aln_query.seq[j] == "-":
                query_gaps = query_gaps + 1
            query_pos = (1 if self.aln_query.strand == 1 else -1) * (j - query_gaps) + q_start

            if self.aln_hit.seq[j] == "-":
                hit_gaps = hit_gaps + 1
            hit_pos = (1 if self.aln_hit.strand == 1 else -1) * (j - hit_gaps) + h_start
            # assert self.aln_query.seq[j] == self.aln_query.seq.replace("-","")[query_pos]
            assert (self.aln_hit.seq[j] == "-") or (
                    self.aln_hit.seq[j] == self.aln_hit.seq.replace("-", "")[hit_pos - h_start])
            aln_pos = AlnPos(j, query_pos, self.aln_query.seq[j], hit_pos, self.aln_hit.seq[j], query_gaps, hit_gaps)
            self.positions.append(aln_pos)

    def aln_len(self):
        return len(self.aln_query)

    def update_hit_coverage(self, hit_len):
        self.hit_coverage = self.aln_len / hit_len

    def update_query_coverage(self, query_len):
        self.query_coverage = self.aln_len / query_len

    def update_identity(self, fn_identity=None):
        if not fn_identity:
            fn_identity = self._fn_identity

        self.identity = self._fn_identity(self.aln_query.seq, self.aln_hit.seq)

    def _fn_identity(self, query_seq, hit_seq):
        aln_length = len(query_seq)
        ident = 0
        for i in range(aln_length):
            if query_seq[i] == hit_seq[i]:
                ident = ident + 1
        if aln_length:
            return float("{0:.2f}".format(1.0 * ident / aln_length))
        else:
            return 0

    def is_gap(self, position):
        return self.aln_hit.is_gap(position) & self.aln_query.is_gap(position)

    def map_pos_hit_query(self, pos_hit):
        pos = [x._q_pos for x in self.positions if x._h_pos == pos_hit and x.qaa != "-"]
        if pos:
            return pos[0]
        else:
            None

    def map_pos_query_hit(self, pos_query):
        pos = [x._h_pos for x in self.positions if x._q_pos == pos_query and x.haa != "-"]
        if pos:
            return pos[0]
        else:
            None

    def map_q_alignment(self):
        return map_alignment(self.aln_query.seq, self.aln_hit.seq)

    def map_h_alignment(self):
        return map_alignment(self.aln_hit.seq, self.aln_query.seq)

    def query_record(self):
        return SeqRecord(Seq(self.aln_query.seq), id=self.aln_query.name)

    def hit_record(self):
        return SeqRecord(Seq(self.aln_hit.seq), id=self.aln_hit.name)

    def __getitem__(self, key):
        _key = key
        if key == -1:
            _key = len(self.positions) - 1
        return [x for x in self.positions if x.pos == _key][0]

    def __repr__(self):
        return "SimpleAlignment(" + self.aln_query.name + " - " + self.aln_hit.name + ")"

    def __str__(self):
        return "\n".join([self.aln_query.name, self.aln_query.seq, self.aln_hit.name, self.aln_hit.seq])
