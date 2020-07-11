import Bio.SeqIO as bpio


class CDHit:
    """
    cd-hit -i proteins.faa -o proteins90.faa -c 0.9 -aL 0.9 -aS 0.9 -g 1 -p 1
    """
    @staticmethod
    def get_clusters(clstr_h_or_file):
        h = clstr_h_or_file if hasattr(clstr_h_or_file, "write") else open(clstr_h_or_file)
        try:
            seq_cluster = {}
            cluster_seqs = {}
            clusters_str = h.read().split(">Cluster")[1:]
            for cluster in clusters_str:                
                clines = [x.strip() for x in cluster.strip().split("\n")]
                cid = int(clines[0])
                seqs = [x.split(">")[1].split("...")[0] for x in clines[1:]]
                cluster_seqs[cid] = seqs
                for s in seqs:
                    seq_cluster[s] = cid
            return seq_cluster,cluster_seqs

        finally:
            h.close()

    @staticmethod
    def unip_cluster(clstr_h_or_file, in_fasta_h_or_file, out_fasta_h_or_file):
        """
        
        :return: 
        """
        h_in = in_fasta_h_or_file if hasattr(in_fasta_h_or_file, "write") else open(in_fasta_h_or_file)
        try:
            h_out = out_fasta_h_or_file if hasattr(out_fasta_h_or_file, "write") else open(out_fasta_h_or_file, "w")
            seq_cluster,cluster_seqs = CDHit.get_clusters(clstr_h_or_file)
            for seqrecord in bpio.parse(h_in, "fasta"):
                cid = seq_cluster[seqrecord.id[:19]]
                seqs = cluster_seqs[cid]
                seq_id = seqrecord.id
                if not seqrecord.id.startswith("sp"):
                    seq_id = sorted(seqs, key=lambda x: 1 if x.startswith("sp") else 0)[-1]
                
                gene = ""
                unip_id = ""
                if "GN=" in seqrecord.description:
                    gene = seqrecord.description.split("GN=")[1].split()[0]
                    unip_id = seqrecord.id
                if seq_id.count("|") == 2:                    
                    unip_id = seq_id.split("|")[1]
                assert "|" not in unip_id,     [seq_id,unip_id , seq_id.count("|")]
                if unip_id or gene:
                    seq_id = gene + "__" + unip_id
                seqrecord.id = seq_id
                seqrecord.name = ""
                assert "|" not in seq_id,     seq_id
                bpio.write(seqrecord, h_out, "fasta")

        finally:
            h_in.close()
            h_out.close()
