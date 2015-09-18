__author__ = 'robswift'

class Writer(object):
    """

    """
    def __init__(self, query):
        self.query = query
        self.out = []

    def write_query(self):
        # query ID, status and reason
        if query.triage:
            query_status = 'Discarded'
        else:
            query_status = 'Retained'
            out.append('Query ID: {id}|Chain count: {chain_count}|Dock count: {dock_count}'.format(id=query.pdb_id,
                chain_count=query.chain_count, dock_count=query.dock_count))
            out.append('Status: {status}'.format(status=query_status))
        if query_status == 'Discarded':
            for reason in query.reasons_to_triage:
            out.append('Reason: {reason}'.format(reason=reason))
        # query chain count and chain lengths
        for seq in query.sequences:
            out.append('Chain_{id}: Length = {length}'.format(id=seq.id, length=len(seq)))
        # query dock count and resname
        for lig in query.dock:
            out.append('Ligand: {res}'.format(res=lig.resname))
        out.append('')

    def write_hit(self):
        pass

    def write_chain(self):
        pass

