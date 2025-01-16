from Bio import Entrez #type:ignore
import time
import os

def collect_ncbi_data(search_term, output_file, email, database, file_format, retmax):
    """
    Query NCBI GenBank database with the given search terms. Collect each entry's ID
    and then downloads queried entries in batches.

    :param search_term: Term(s) used to query the database.
    :param output_file: Where GenBank data will be stored.
    :param email: Required by NCBI
    :param database: Selected NCBI database to query.
    :param file_format: Selected data format.
    :param retmax: Maximum number of entries to be returned.

    :output metagenome_data.gb:
    """

    Entrez.email = email

    try:
        stream = Entrez.esearch(db=database, term=search_term, retmax=retmax, usehistory="y", idtype="acc")
        search_results = Entrez.read(stream)
        stream.close()

        ids = search_results.get("IdList", [])

        with open(output_file, "w") as file:
            """
            Without an API key, NCBI limits users to 3 queries per second.
            Thus, search results are broken into batches and downloaded batchwise.
            """

            for i in range(0, len(ids), 3):
                batch_ids = ids[i: i+3]
                try:
                    batch_stream = Entrez.efetch(db=database, id=",".join(batch_ids), rettype=file_format, retmode="text")
                    file.write(batch_stream.read())
                    batch_stream.close()
                    # Pause 1 second to respect rate limit
                    time.sleep(1)
                except Exception as e:
                    print(f"Error: {e}")

    except Exception as e:
        print(f"Error: {e}")

def run_collector():

    cd = os.path.dirname(os.path.abspath(__file__))
    output_file = os.path.join(cd, '../../data/metagenome_data.gb')

    search_term = "metagenome[ORGN] AND contig"
    email = "jg198@student.london.ac.uk"
    database = "nucleotide"
    file_format = "gb"
    # Total number of results for current search term, as of January 2025
    retmax = 659
    collect_ncbi_data(search_term, output_file, email, database, file_format, retmax)

"""
To run the data collector:
1. Uncomment the function call below
2. Open a Terminal instance at the root directory
3. Ensure you are working with the correct interpreter
4. In the terminal, type "python3 -m src.data_pipeline.data_collector"
5. Wait for the process to complete or end it early by pressing "Ctrl+C"
"""
#run_collector()