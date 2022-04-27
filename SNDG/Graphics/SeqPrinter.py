from dna_features_viewer import GraphicFeature

from SNDG.Graphics.CustomBiopythonTranslator import CustomBiopythonTranslator


class SeqPrinter:

    def __init__(self, seq_record, start=0, end=None, translator=None):
        self.seq_record = seq_record[:]
        self._start = start
        self._end = end
        self.translator = translator if translator else CustomBiopythonTranslator.buildDefault()
        self.gfeatures = []

    def start(self):
        return self._start if self._start >=0 else 0

    def end(self):
        if self._end:
            end = self._end if self._end < len(self.seq_record) else len(self.seq_record)
        else:
            end = len(self.seq_record)
        return end

    def print(self, ax=None, with_ruler=True, gfeatures=None, add_contig=True,
              figure_width=15,file_path=None):

        # mct.contig = contig_record.id
        record = self.translator.translate_record(self.seq_record)
        if add_contig:
            r = GraphicFeature(
                start=0, end=len(self.seq_record), strand=+1, color="yellow", label=self.seq_record.id)
            record.features.append(r)

        if gfeatures:
            for gf in gfeatures:
                record.features.append(gf)

        cropped_record = record.crop((self.start(), self.end()))

        if not ax:
            ax, _ = cropped_record.plot(figure_width=figure_width, with_ruler=with_ruler)
        else:
            cropped_record.plot(ax=ax, figure_width=figure_width, with_ruler=with_ruler)
        if file_path:
            ax.figure.savefig(file_path)
