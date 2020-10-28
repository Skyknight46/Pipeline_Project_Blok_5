def fasta_adder(originele_fasta, nieuwe_fasta, aanpasbare_file_location):
    """
    Voegt de sequenties van het naar fasta geconverteerde hmmsearch
    stockholm bestand toe aan het originele fasta bestand
    :param originele_fasta: path naar de ClustalOmega inputfasta van
    een bepaalde iteratie
    :param nieuwe_fasta: path naar de hmmSearch output van dezelfde iteratie
    :param aanpasbare_file_location: Zelf gespecificeerde locatie waar de
    samengevoegde sequenties in komen te staan
    :return: pad naar file met samengevoegde sequenties (inputsequenties bij
     bepaalde iteraties)
    """
    originele_accessie_set = set()
    nieuwe_accessie_set = set()
    nieuwe_dict = {}
    # loopt over regels van FASTA bestand en voegt de
    # accessiecodes toe aan een set
    with open(originele_fasta, 'r') as open_file:
        for line in open_file:
            if line.startswith(">"):
                originele_accessie_set.add(line.strip()
                                           .replace(">", "").split(" ")[0])
    # loopt over regels van een ander FASTA bestand,, voegt a
    # ccessiecodes toe aan een set,
    # splitst headers en sequenties en voegt deze toe aan een dictionary
    with open(nieuwe_fasta, 'r') as open_file:
        for line in open_file:
            if line.startswith(">"):
                accessie = line.strip().replace(">", "").split("/")[0]
                key = line.strip()
                nieuwe_accessie_set.add(accessie)
                nieuwe_dict[key] = ""
            else:
                nieuwe_dict[key] += line.strip().upper()
    # maakt een nieuwe set met de headers die wel in het het 2de FASTA
    # bestand zitten maar niet in het 1ste
    verschil = nieuwe_accessie_set.difference(originele_accessie_set)
    # voegt de headers en sequenties die bij de accessiecodes in
    # 'verschil' horen toe bij een FASTA bestand
    returnfile = aanpasbare_file_location
    with open(aanpasbare_file_location, "a") as open_file:
        for accessiecode in list(verschil):
            for key, value in nieuwe_dict.items():
                if accessiecode in key:
                    open_file.write(key)
                    open_file.write("\n")
                    open_file.write(nieuwe_dict[key])
                    open_file.write("\n")
        return returnfile
