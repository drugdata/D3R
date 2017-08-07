__author__ = 'robswift'

import os

class Writer(object):
    """
    Base writer class
    """
    def __init__(self, out_dir):
        self.query = None
        self.out_dir = out_dir
        self.handle = None

    def write_query_header(self):
        out = []
        out.append('Query ID: {id}|Sequence count: {seq_count}|Dock count: {dock_count}'.format(
            id=self.query.pdb_id, seq_count=self.query.sequence_count, dock_count=self.query.dock_count))
        for l in out: self.handle.write("%s\n" % l)

    def write_query_status(self):
        """
        Check set and write the status of a query.
        """
        out = []
        if self.query.triage:
            query_status = 'Discarded'
        else:
            query_status = 'Retained'
        out.append('Status: {status}'.format(status=query_status))
        if query_status == 'Discarded':
            for reason in self.query.reasons_to_triage:
                out.append('Reason: {reason}'.format(reason=reason))
        for l in out: self.handle.write("%s\n" % l)

    def write_query_chains(self):
        out = []
        for seq in self.query.sequences:
            out.append('Chain_{id}: Length = {length}'.format(id=seq.id.split('_')[1], length=len(seq)))
        for l in out: self.handle.write("%s\n" % l)

    def write_query_ligands(self):
        if self.query.dock_count == 0:
            pass
        else:
            out = []
            for lig in self.query.dock:
                out.append('Ligand: {res}'.format(res=lig.resname))
            for l in out: self.handle.write("%s\n" % l)

    def write_query_ph(self):
        out = []
        if self.query.exp_ph:
            out.append('Crystallization pH: {ph}'.format(ph=self.query.exp_ph))
        else:
            out.append('Crystallization pH: N/A')
        for l in out: self.handle.write("%s\n" % l)

    def write_hit_header(self, hit):
        out = []
        out.append('  Hit ID: {id}|Sequence count: {seq_count}|Chain count: {chain_count}|Dock count:'
                   ' {dock_count}'.format(id=hit.pdb_id, seq_count=hit.sequence_count, chain_count=hit.chain_count,
                                          dock_count=hit.dock_count))
        for l in out: self.handle.write("%s\n" % l)

    def write_hit_status(self, hit):
        out = []
        if hit.triage:
            hit_status = 'Discarded'
        elif hit.retain:
            hit_status = 'Elected for docking'
        else:
            hit_status = 'Docking candidate'
        out.append('  Status: {stat}'.format(stat=hit_status))
        if hit_status == 'Discarded':
            for reason in hit.reasons_to_triage:
                out.append('  Reason: {reason}'.format(reason=reason))
        elif hit_status == 'Elected for docking':
            for reason in hit.reasons_to_retain:
                out.append('  Reason: {reason}'.format(reason=reason))
        for l in out: self.handle.write("%s\n" % l)

    def write_hit_chains(self, hit):
        out = []
        hit.sequences.sort(key=lambda seq: seq.blast_hit, reverse=True)
        for seq in hit.sequences:
            if seq.blast_hit:
                for qa in seq.query_alignments:
                    out.append('  Chain_{id}_{no}:  Length = {length}|%identity = {identity:.2} '
                               '(with Chain_{query_id})|%coverage = {coverage:.2} '
                               '(with Chain_{query_id})'.format(id=seq.hit_chain_id, no=seq.hit_sequence_id,
                                                                length=len(seq.seq_record), identity=qa.identity,
                                                                coverage=qa.coverage, query_id=qa.query_chain_id))
            else:
                out.append('  Chain_{id}_{no}: Length = {length}|Not a BLAST hit'.format(id=seq.hit_chain_id,
                                                                                          no=seq.hit_sequence_id,
                                                                                          length=len(seq.seq_record)))

        for l in out: self.handle.write("%s\n" % l)

    def write_hit_ligands(self, hit):
        if hit.dock_count == 0:
            pass
        else:
            out = []
            for lig in hit.dock:
                line = []
                line.append('  Ligand: {res}'.format(res=lig.resname))
                for mcss in lig.mcsss:
                    line.append('|MCSS size: {s} (with {l})'.format(s=mcss.size, l=mcss.reference))
                out.append(''.join(line))
            for l in out: self.handle.write("%s\n" % l)

    def write_hit_method(self, hit):
        out = []
        out.append('  Method: {method}|Resolution: {res}'.format(method=hit.exp_method, res=hit.resolution))
        for l in out: self.handle.write("%s\n" % l)

    def write_blank_lines(self, number):
        out = []
        for i in range(number):
            out.append('')
        for l in out: self.handle.write("%s\n" % l)

    def write_hit(self, hit):
        self.write_hit_header(hit)
        self.write_hit_status(hit)
        self.write_hit_chains(hit)
        self.write_hit_ligands(hit)
        self.write_hit_method(hit)
        self.write_blank_lines(1)


class WriteText(object):
    """
    This is the new txt writer that meets the format that Shuai wanted.
    """
    def __init__(self, out_dir):
        self.query = None
        self.out_dir = out_dir
        self.handle = None
        self.most_similar = []
        self.least_similar = []
        self.highest_res = []
        self.apo = []
        self.highest_tanimoto = []
        self.reasons = {
            1: 'The BLAST hit is bound to the ligand with the largest maximum common substructure',
            2: 'The BLAST hit is bound to the ligand with the smallest maximum common substructure',
            3: 'The BLAST hit is the highest resolution holo structure',
            4: 'The BLAST hit is the highest resolution apo structure',
            5: 'The BLAST hit is bound to the ligand with the highest tanimoto score in structural similarity',
        }

    def write_txt(self, query):
        self.reinitialize(query)
        self.categorize()
        self.open_file()
        self.write_query()
        self.write_hits()
        self.close_file()

    def write_query(self):
        self.write_query_header()
        self.write_query_ph()
        self.write_query_ligands()

    def write_query_header(self):
        out = []
        out.append('query, {id}'.format(id=self.query.pdb_id))
        for l in out: self.handle.write("%s\n" % l)

    def write_query_ph(self):
        out=[]
        if self.query.exp_ph:
            out.append('ph, {ph}'.format(ph=self.query.exp_ph))
        else:
            out.append('ph, N/A')
        for l in out: self.handle.write("%s\n" % l)

    def write_query_ligands(self):
        if self.query.dock_count == 0:
            pass
        else:
            out = []
            for lig in self.query.dock:
                out.append('ligand, {res}'.format(res=lig.resname))
                out.append('inchi, {inchi}'.format(inchi=lig.inchi))
                #modified by sliu 08/08 
                out.append('size, {size}'.format(size=lig.size))
                out.append('rotatable_bond, {rot}'.format(rot=lig.rot))
                #out.append('heavy size, {heavy_size}'.format(heavy_size=lig.heavy))
            for l in out: self.handle.write("%s\n" % l)

    def write_hits(self):
        #here change to only write out the top 10 largest
        for hit in self.most_similar[:10]:
            self.write_largest(hit)
        for hit in self.least_similar:
            self.write_smallest(hit)
        #sliu 01092017 change to top ten in holo structures
        for hit in self.highest_res[:10]:
            self.write_holo(hit)
        for hit in self.apo:
            self.write_apo(hit)
        for hit in self.highest_tanimoto[:10]:
            self.write_highest_tanimoto(hit)

    def write_largest(self, hit):
        out = []
        #only write out the dockable mcss in the first ordered chain
        largest_lig = hit.dock[hit.largest_index[0]]
        for mcss in largest_lig.mcsss:
            out.append("LMCSS, {pdb_id}, {lig}, chain: {chain_id}, (size: {lig_size}, mcss_size: {mcss_size}, resolution: {resolution}) ".format(pdb_id=hit.pdb_id, lig=mcss.test, chain_id=hit.largest_mcss_chain[0], lig_size= largest_lig.size, mcss_size= mcss.size, resolution= "%4.4s"%hit.resolution ))
        for l in out: self.handle.write("%s\n" % l)

    def write_smallest(self, hit):
        out = []
        smallest_lig = hit.dock[hit.smallest_index[0]]
        for mcss in hit.dock[hit.smallest_index[0]].mcsss:
            #out.append("SMCSS, {pdb_id}, {lig} ".format(pdb_id=hit.pdb_id, lig=mcss.test ))
            out.append("SMCSS, {pdb_id}, {lig}, chain: {chain_id}, (size: {lig_size}, mcss_size: {mcss_size}, resolution: {resolution}) ".format(pdb_id=hit.pdb_id, lig=mcss.test, chain_id=hit.smallest_mcss_chain[0], lig_size= smallest_lig.size, mcss_size= mcss.size, resolution= "%4.4s"%hit.resolution ))
        for l in out: self.handle.write("%s\n" % l)

    def write_holo(self, hit):
        out = []
        largest_lig = hit.dock[hit.largest_index[0]] 
        for mcss in largest_lig.mcsss:
            out.append("hiResHolo, {pdb_id}, {lig}, chain: {chain_id}, (resolution: {resolution})".format(pdb_id=hit.pdb_id, lig=mcss.test, chain_id=hit.largest_mcss_chain[0], resolution= "%4.4s"%hit.resolution))
#        for lig in hit.dock[:1]:
#            for mcss in lig.mcsss:
#                out.append("hiResHolo, {pdb_id}, {lig}".format(pdb_id=hit.pdb_id, lig=mcss.test))
        for l in out: self.handle.write("%s\n" % l)

    def write_apo(self, hit):
        out = []
        out.append("hiResApo, {pdb_id}".format(pdb_id=hit.pdb_id))
        for l in out: self.handle.write("%s\n" %l)
    def write_highest_tanimoto(self, hit):
        out = []
        #only write out the dockable mcss in the first ordered chain
        largest_lig = hit.dock[hit.highest_tanimoto_index[0]]
        for mcss in largest_lig.mcsss:
            out.append("hiTanimoto, {pdb_id}, {lig}, chain: {chain_id}, (tanimoto_similarity: {tanimoto}, resolution: {resolution})".format(pdb_id=hit.pdb_id, lig=mcss.test, chain_id=hit.highest_tanimoto_chain[0], tanimoto= "%4.4s"%(mcss.tanimoto), resolution= "%4.4s"%hit.resolution  ))
        for l in out: self.handle.write("%s\n" % l)

    def reinitialize(self, query):
        self.most_similar = []
        self.least_similar = []
        self.highest_res = []
        self.apo = []
        self.query = query

    def open_file(self):
        f = os.path.join(self.out_dir, '{stem}.txt'.format(stem=self.query.pdb_id))
        self.handle = open(f, 'w')

    def close_file(self):
        self.handle.close()

    def categorize(self):
        for hit in self.query.hits:
            if hit.retain:
                for reason in hit.reasons_to_retain:
                    if reason == self.reasons[1] and hit.pdb_id:
                        self.most_similar.append(hit)
                    if reason == self.reasons[2] and hit.pdb_id:
                        self.least_similar.append(hit)
                    if reason == self.reasons[3] and hit.pdb_id:
                        self.highest_res.append(hit)
                    if reason == self.reasons[4] and hit.pdb_id:
                        self.apo.append(hit)
                    #add tanimoto 09/01 sliu
                    if reason == self.reasons[5] and hit.pdb_id:
                        self.highest_tanimoto.append(hit)
        #here rank the list by the resolution
        self.most_similar.sort(key=lambda hit:(float(hit.resolution)))
        self.least_similar.sort(key=lambda hit:(float(hit.resolution)))
        self.highest_res.sort(key=lambda hit:(float(hit.resolution)))
        self.apo.sort(key=lambda hit:(float(hit.resolution)))
        self.highest_tanimoto.sort(key=lambda hit:(float(hit.resolution)))


# class WriteTxt(Writer):
#    """
#    Writes a txt file that describes a query and the corresponding blast hits elected for docking. This is the old
#    class that is less useful for Shuai.
#    """
#    def __init__(self, *args, **kwargs):
#        super(WriteTxt, self).__init__(*args, **kwargs)
#        self.most_similar = []
#        self.least_similar = []
#        self.highest_res = []
#        self.apo = []
#        self.reasons = {
#            1: 'The BLAST hit is bound to the ligand with the largest maximum common substructure',
#            2: 'The BLAST hit is bound to the ligand with the smallest maximum common substructure',
#            3: 'The BLAST hit is the highest resolution structure',
#            4: 'The BLAST hit is the highest resolution apo structure',
#        }
#
#    def write_txt(self, query):
#        self.reinitialize(query)
#        self.categorize()
#        self.open_file()
#        self.write_query()
#        self.write_hits()
#        self.close_file()
#
#    def write_query(self):
#        self.write_query_header()
#        self.write_query_chains()
#        self.write_query_ligands()
#        self.write_blank_lines(1)
#
#    def write_hits(self):
#        for hit in self.most_similar:
#            self.write_hit(hit)
#        for hit in self.least_similar:
#            self.write_hit(hit)
#        for hit in self.highest_res:
#            self.write_hit(hit)
#        for hit in self.apo:
#            self.write_hit(hit)
#
#    def reinitialize(self, query):
#        self.most_similar = []
#        self.least_similar = []
#        self.highest_res = []
#        self.apo = []
#        self.query = query
#
#    def open_file(self):
#        f = os.path.join(self.out_dir, '{stem}.txt'.format(stem=self.query.pdb_id))
#        self.handle = open(f, 'w')
#
#    def close_file(self):
#        self.handle.close()
#
#    def categorize(self):
#        assigned = []
#        for hit in self.query.hits:
#            if hit.retain:
#                for reason in hit.reasons_to_retain:
#                    if reason == self.reasons[1] and hit.pdb_id not in assigned:
#                        self.most_similar.append(hit)
#                        assigned.append(hit.pdb_id)
#                    if reason == self.reasons[2] and hit.pdb_id not in assigned:
#                        self.least_similar.append(hit)
#                        assigned.append(hit.pdb_id)
#                    if reason == self.reasons[3] and hit.pdb_id not in assigned:
#                        self.highest_res.append(hit)
#                        assigned.append(hit.pdb_id)
#                    if reason == self.reasons[4] and hit.pdb_id not in assigned:
#                        self.apo.append(hit)
#                        assigned.append(hit.pdb_id)

class WriteLog(Writer):
    """
    Writes a log file that describes a query and the corresponding blst hits
    """
    def __init__(self, *args, **kwargs):
        super(WriteLog, self).__init__(*args, **kwargs)
        self.open_file()

    def open_file(self):
        f = os.path.join(self.out_dir, 'blastnfilter.log')
        self.handle = open(f, 'a')

    def close_file(self):
        self.handle.close()

    def write_log(self, query):
        self.query = query
        self.write_query()
        self.write_hits()

    def write_query(self):
        self.write_query_header()
        self.write_query_status()
        self.write_query_chains()
        self.write_query_ligands()
        self.write_query_ph()
        self.write_blank_lines(1)

    def write_hits(self):
        for hit in self.query.hits:
            self.write_hit(hit)
