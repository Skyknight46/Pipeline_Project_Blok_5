from Bio import Entrez, GenBank
import ssl
import time

# In de entrez_id_search en entrez_lineage functie, vervang PLACEHOLDER door je eigen 
# mailadres zodat je een id krijgt bij de NCBI



def entrez_id_search(accessiecode):  # mogen maar 3 request per seconde
    """
    Zoekt in de protein db van ncbi doormiddel van entrez naar het id
    zodat dat gebruikt kan worden om het genbank bestand op te zoeken
    en de taxonomy te gebruiken
    :param accessiecode: Accessiecode van huidige entry
    :return result: Dictionary met als key de accessiecode en als value de lineage
    """
    # print(accessiecode, "in entrez_id_search")
    ssl._create_default_https_context = ssl._create_unverified_context
    time.sleep(0.5)  # voorkomt overload op ncbi server, anders ban :)
    Entrez.email = "PLACEHOLDER"  
    handle = Entrez.esearch(db="protein", term=accessiecode)
    # zoekt in de protein db naar de accessiecode
    record = Entrez.read(handle)
    id_access = record["IdList"]
    result = entrez_lineage(id_access)
    return result  # wordt doorgegeven aan functie taxonomy in gap_results.py


def entrez_lineage(id_access):
    """
    Zoekt in de protein db met het id naar het genbank bestand
    waar de taxonomy in staat.
    :param id_access: unieke ID voor elke accessiecode
    :return: Dictionary met als key de accessiecode en als value de lineage
    """
    try:
        Entrez.email = "PLACEHOLDER"  # verplicht
        filename = "genbankfile.txt"
        with Entrez.efetch(db="protein", id=id_access, rettype="gb",
                           retmode="text") as net_handle:
            with open(filename, "w") as output_handle:
                output_handle.write(net_handle.read())
                # print("Genbank file opgeslagen")
        regions = {}
        with open(filename, "r") as handle:
            region = False
            record = GenBank.read(handle)
            count = 0
            count1 = 0
            key = record.locus
            for i in record.features:
                if "Region" in str(record.features[count]) and count1 != 0:
                    regions[key] += str(record.features[count]).\
                        replace("\n","").replace(
                        " " * 20, " ")
                if "Region" in str(record.features[count]) and count1 == 0:
                    region = True   # kijkt of er regions in het gebank
                    # bestand zitten
                    regions[key] = str(record.taxonomy) + " " + \
                        str(record.features[count]).replace("\n", "")\
                        .replace(" " * 20, " ")
                    count1 += 1
                count += 1
            if not region:
                regions[key] = str(record.taxonomy)
            return regions
    except ValueError:  # als er geen taxonomy bekend is
        return "No taxonomy known"
    except KeyError:
        return "No key"

# geholpen websites:
# https://biopython.org/docs/dev/api/Bio.GenBank.html
# https://biopython-tutorial.readthedocs.io/en/latest/notebooks/
# 09%20-%20Accessing%20NCBIs%20Entrez%20databases.html
# https://biopython.org/docs/1.75/api/Bio.GenBank.Record.html lijst
# met functies
