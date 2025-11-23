README ASSIGNMENT1 ASB

BUSCAR SEQUÊNCIAS DO GENBANK:
db="nucleotide" #Define a base de dados do NCBI
accessions="JQ307171,JQ307172,JQ307173,JQ307174,JQ307175,JQ307176,JQ307177,JQ307178,JQ307179,JQ307180,JQ307181,JQ307182,JQ307183,JQ307184,JQ307185,JQ307186,JQ307187,JQ307188,JQ307189,JQ307190,JQ307191,JQ307192,JQ307193,JQ307194,JQ307195,JQ307196,JQ307197,JQ307198,JQ307199,JQ307200"
#Lista de Acession Numbers correspondentes ao gene COI
wget -qO COI_sequences.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$db&id=$accessions&rettype=fasta&retmode=text"
#Retorna as sequências de DNA em formato FASTA e grava no ficheiro COI_sequences.fasta
#Fizemos o mesmo para os genes white e cad com os seus respetivos accession numbers

ALINHAMENTO DAS SEQUÊNCIAS
#Usamos o mafft para alinhar as sequências e acompanhamos o resultado do alinhamento com o aliview
mafft COI_sequences.fasta > align_COI.fasta
mafft CAD_sequences.fasta > align_CAD.fasta
mafft white_sequences.fasta > align_white.fasta

ALteramos os nomes das amostras nos diferentes genes para que ficassem iguais tornando a concatenação possível


FILTRAR FASTA (SCRIPT PYTHON) #Este script elimina os Accession Numbers presentes nos nomes das Amostras


from Bio import SeqIO
import re

def clean_headers(input_file, output_file):
    """
    Mantém apenas a parte do cabeçalho que começa com 'Anopheles'
    """
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            # Pega só o texto a partir de 'Anopheles'
            new_header = re.search(r"Anopheles.*", record.description)
            if new_header:
                record.id = new_header.group().split()[0]
                record.description = new_header.group()
            SeqIO.write(record, out_handle, "fasta")

# Executa se for chamado diretamente
if __name__ == "__main__":
    input_fasta = "C:/Users/rgs12/Downloads/align_white2.fasta"
    output_fasta = "C:/Users/rgs12/Downloads/filtered_white.fasta"
    
    clean_headers(input_fasta, output_fasta)
    print(f"Novo ficheiro criado: {output_fasta}")
#Foi feito o mesmo para o resto dos genes


AGRUPAR NOMES (SCRIPT PYTHON)
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def process_fasta(input_file, output_file):
    replacements = {
        # An. galvaoi
        r'Anopheles galvaoi voucher PR19_2_101.*': 'An. galvaoi',
        r'Anopheles galvaoi voucher SP18_111.*': 'An. galvaoi',
        r'Anopheles galvaoi voucher SP66_20_1.*': 'An. galvaoi',
with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, 'fasta'):
            new_id = record.id
            new_description = record.description
            
            # Processa tanto o ID quanto a descrição
            for pattern, replacement in replacements.items():
                if re.match(pattern, new_id):
                    new_id = replacement
                    break
                    
            for pattern, replacement in replacements.items():
                if re.match(pattern, new_description):
                    new_description = replacement
                    break
            
            # Cria um novo registro com os dados modificados
            new_record = SeqRecord(
                record.seq,
                id=new_id,
                description=new_description,
                name=''  # Limpa o nome para evitar duplicações
            )
            
            # Escreve o registro modificado no arquivo de saída
            SeqIO.write(new_record, out_handle, 'fasta')



process_fasta('C:/Users/rgs12/Downloads/filtered_coi.fasta', 'coi_agrupado.fasta')


CONCATENAR FASTA 
#Usamos o FASconCAT-G-master e movemos os 3 fasta para este diretório para concatenar os ficheiros 
wget https://github.com/PatrickKueck/FASconCAT-G/archive/master.zip 
unzip master.zip   
cd FASconCAT-G-master 
perl FASconCAT-G_v1.06.1.pl 
#Este software gerou um arquivo fasta já concatenado

CONVERTER FASTA EM NEXUS
#Usamos o mega para fazer a conversão do ficheiro e gerou o ficheiro: 

Árvore BAYESIANA
mb concatenado.nexus 

ANÁLISE DA ÁRVORE
figtree

GITHUB
git init
git remote add origin git@github.com:rodrigosilva1936/ASB.git
git add .
git commit -m "Adiciona scripts iniciais e dados"
git branch -M main
git push -u origin main
hash commit : 0c014ad1ff13c865e078a3bf32a010df0850931c


Resumo de Softwares e Aplicações Usadas
mafft (alinhamento)
aliview (análise das sequências alinhadas)
FASconCAT-G-master (concatenação)
Mega (conversão para nexus)
MrBayes (criar árvore bayesiana)
figtree (Análise da Árvore Filogética)