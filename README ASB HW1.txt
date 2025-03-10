Códigos usados:

nano fetch_fasta.sh - usamos este comando para abrir o ficheiro fetch_fasta.sh e colámos o código que usámos para a realização deste homework.

Bloco de Comandos gravado posteriormente dentro do comando nano:
db="nucleotide"                                                            
query="Psammodromus algirus[organism], cytb[gene]"
search_result=$(wget -qO - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=$db&term=$query&usehistory=y")
webenv=$(echo "$search_result" | grep -oP '(?<=<WebEnv>).*?(?=</WebEnv>)')
query_key=$(echo "$search_result" | grep -oP '(?<=<QueryKey>).*?(?=</QueryKey>)')
wget -qO - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$db&WebEnv=$webenv&query_key=$query_key&rettype=fasta"


db="nucleotide" - Define a Base de Dados a ser consultado (nucleotide que contém sequências de DNA)
query="Psammodromus algirus[organism], cytb[gene]" - Define a pesquisa (sequências do gene cytb no organismo Psammodromus algirus

search_result=$(wget -qO - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=$db&term=$query&usehistory=y") 
- wget -qO - "URL" faz uma requisição HTTP silenciosa e retorna o resultado no terminal, db=$db (Banco de dados nucleotide), term=$query (pesquisa usando o nome do organismo e gene), usehistory = y - Habilita o History Server, permitindo recuperar os resultados futuramente sem precisar de procurar outravez. O resultado é armazenado na variável search_result

webenv=$(echo "$search_result" | grep -oP '(?<=<WebEnv>).*?(?=</WebEnv>)') - Extrai o valor dentro da tag <WebEnv>...</WebEnv>
query_key=$(echo "$search_result" | grep -oP '(?<=<QueryKey>).*?(?=</QueryKey>)') - Extrai o valor dentro da tag <QueryKey>...</QueryKey>

wget -qO - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$db&WebEnv=$webenv&query_key=$query_key&rettype=fasta" - efetch.fcgi é o serviço da API usado para recuperar dados de uma pesquisa anterior, db=$db (Banco de dados nucleotide),  WebEnv=$webenv - usa o identificador de sessão recuperado anteriormente, query_key=$query_key  - Usa o identificador da pesquisa específica, rettype=fasta → Define o formato de saída como FASTA.

Após o código:
chmod +x fetch_fasta.sh - Torna o ficheiro executável
./fetch_fasta.sh > resultado.fasta - converte o ficheiro para um ficheiro fasta

Instruções para o commit - 
git add resultado.fasta
git commit -m "Adicionei o Homework"
git remote add origin
git push origin main 
git log --oneline (A key é 18fc1e4)

Trabalho realizado por:
Rodrigo Silva 202200138
Núria Machado 202300237
Beatriz Venâncio 202300392
Erica Martins 202300387





















