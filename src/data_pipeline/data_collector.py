from Bio import Entrez#type:ignore

class DataCollector:
    def __init__(self, email, output_dir="data", batch_size=10, db="nucleotide"):
        self.db = db
        self.output_dir = output_dir
        self.batch_size = batch_size
        Entrez.email = email


    def search_metagenomes(self, search_term, retmax=100):
        try:
            handle = Entrez.esearch(db=self.db, term=search_term, retmax=retmax)
            record = Entrez.read(handle)
            return record.get("IdList", [])
        except Exception as e:
            print(f"error: {e}")
            return []

    
    def fetch_metagenomes(self, ids, file_format="gb"):
        # TODO - do
        pass

    def save_metagenomes(self):
        # TODO - do
        pass

        
if __name__ == "__main__":
    # Data Collector Parameters
    email = "jg198@student.london.ac.uk"
    search_term = "metagenome[ORGN] AND contig"
    output_dir = "metagenome_data" # TODO - Align file
    database = "nucleotide"
    collector = DataCollector(email,
                              output_dir=output_dir,
                              search_term=search_term,
                              db=database)


