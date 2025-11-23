import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def process_fasta(input_file, output_file):
    replacements = {
        # An. galvaoi
        r'Anopheles galvaoi voucher PR19_2_101.*': 'An. galvaoi',
        r'Anopheles galvaoi voucher SP18_111.*': 'An. galvaoi',
        r'Anopheles galvaoi voucher SP66_20_1.*': 'An. galvaoi',
        
        # An. rangeli
        r'Anopheles rangeli voucher AC15_04.*': 'An. rangeli',
        r'Anopheles rangeli voucher AC18_110.*': 'An. rangeli',
        r'Anopheles rangeli voucher RO18_8_2.*': 'An. rangeli',
        
        # An. konderi A
        r'Anopheles konderi voucher AP15_11.*': 'An. konderi A',
        r'Anopheles konderi voucher AP21_43.*': 'An. konderi A',
        r'Anopheles konderi voucher AP25_11_24.*': 'An. konderi A',
        r'Anopheles konderi voucher AP25_1_100.*': 'An. konderi A',
        r'Anopheles konderi voucher FSP-USP:AP15_11.*': 'An. konderi A',
        r'Anopheles konderi voucher FSP-USP:AP21_43.*': 'An. konderi A',
        r'Anopheles konderi voucher FSP-USP:AP25_11_24.*': 'An. konderi A',
        r'Anopheles konderi voucher FSP-USP:AP25_1_100.*': 'An. konderi A',
        
        # An. konderi s. s.
        r'Anopheles konder voucher PR06_2_13.*': 'An. konderi s. s.',
        r'Anopheles konder voucher PR14_1_9.*': 'An. konderi s. s.',
        r'Anopheles konder voucher PR14_3_17.*': 'An. konderi s. s.',
        r'Anopheles konder voucher PR14_9_108.*': 'An. konderi s. s.',
        r'Anopheles konder voucher RO18_1_6.*': 'An. konderi s. s.',
        r'Anopheles konderi voucher FSP-USP:PR06_2_13.*': 'An. konderi s. s.',
        r'Anopheles konderi voucher FSP-USP:PR14_1_9.*': 'An. konderi s. s.',
        r'Anopheles konderi voucher FSP-USP:PR14_3_17.*': 'An. konderi s. s.',
        r'Anopheles konderi voucher FSP-USP:PR14_9_108.*': 'An. konderi s. s.',
        r'Anopheles konderi voucher FSP-USP:RO18_1_6.*': 'An. konderi s. s.',
        
        # An. konderi B
        r'Anopheles konderi voucher AC18_16.*': 'An. konderi B',
        r'Anopheles konderi voucher FSP-USP:AC18_16.*': 'An. konderi B',
        
        # An. oswaldoi A
        r'Anopheles oswaldoi s. l. AC18 voucher AC18_102.*': 'An. oswaldoi A',
        r'Anopheles oswaldoi s. l. AC18 voucher AC18_107.*': 'An. oswaldoi A',
        r'Anopheles oswaldoi s. l. PA_15 voucher PA_15_C1F2.*': 'An. oswaldoi A',
        r'Anopheles oswaldoi s. l. PA_15 voucher PA_15_C2F4.*': 'An. oswaldoi A',
        r'Anopheles oswaldoi voucher AC18_102.*': 'An. oswaldoi A',
        r'Anopheles oswaldoi AC18_107.*': 'An. oswaldoi A',
        r'Anopheles oswaldoi s. l. PA_15_C1F2.*': 'An. oswaldoi A',
        r'Anopheles oswaldoi s. l. PA_15_C2F4.*': 'An. oswaldoi A',
        r'Anopheles oswaldoi s. l. PA_15_C2F7.*': 'An. oswaldoi A',
        
        # An. oswaldoi ss
        r'Anopheles oswaldoi voucher FSP-USP:ES08_11_07.*': 'An. oswaldoi ss',
        r'Anopheles oswaldoi voucher P03_06SP.*': 'An. oswaldoi ss',
        r'Anopheles oswaldoi voucher SP22_9.*': 'An. oswaldoi ss',
        r'Anopheles oswaldoi voucher ES08_11_07.*': 'An. oswaldoi ss',
        
        # An. dunhami
        r'Anopheles dunhami voucher BRAM13_06.*': 'An. dunhami',
        r'Anopheles dunhami voucher BRAM13_07SP.*': 'An. dunhami',
        r'Anopheles dunhami voucher BRAM13_07SP clone 1.*': 'An. dunhami',
        r'Anopheles dunhami voucher BRAM13_07SP clone 2.*': 'An. dunhami',
        r'Anopheles dunhami voucher BRAM13_113.*': 'An. dunhami',
        
        # An. evansae
        r'Anopheles evansae voucher PR19_10_104.*': 'An. evansae',
        r'Anopheles evansae voucher SP12_28.*': 'An. evansae',
        r'Anopheles evansae voucher SP12_44.*': 'An. evansae',
        r'Anopheles evansae voucher SP18_106.*': 'An. evansae',
        r'Anopheles evansae voucher SP18_27.*': 'An. evansae',
        r'Anopheles evansae voucher VP06_7_4.*': 'An. evansae',
        
        # An. nuneztovari
        r'Anopheles nuneztovari voucher RO1_107.*': 'An. nuneztovari',
        r'Anopheles nuneztovari voucher RO20_02_03.*': 'An. nuneztovari',
        r'Anopheles nuneztovari voucher RO2_13.*': 'An. nuneztovari',
        r'Anopheles nuneztovari voucher RO4_02.*': 'An. nuneztovari',
        
        # An. goeldii
        r'Anopheles goeldii voucher FSP-USP:BRAM03_01.*': 'An. goeldii',
        r'Anopheles goeldii voucher FSP-USP:BRAM22_101.*': 'An. goeldii',
        r'Anopheles goeldii voucher PA7_02_02.*': 'An. goeldii',
        r'Anopheles goeldii voucher PA7_03_08.*': 'An. goeldii',
        r'Anopheles goeldii voucher PA7_04_03.*': 'An. goeldii',
        r'Anopheles goeldii voucher PA7_17_02.*': 'An. goeldii',
        r'Anopheles goeldii voucher BRAM03_01.*': 'An. goeldii',
        r'Anopheles goeldii voucher BRAM22_101.*': 'An. goeldii',
        
        # An. strodei
        r'Anopheles strodei voucher ES09_1.*': 'An. strodei',
        r'Anopheles strodei voucher MG30_102.*': 'An. strodei',
        r'Anopheles strodei voucher SPR04_07 clone 7.*': 'An. strodei',
        r'Anopheles strodei voucher SPR04_07.*': 'An. strodei',
        r'Anopheles strodei voucher VP06_05_01.*': 'An. strodei',
        
        # An. albertoi
        r'Anopheles strodei voucher MG07_12_4.*': 'An. albertoi',
        r'Anopheles strodei voucher MG07_7_10.*': 'An. albertoi',
        r'Anopheles strodei voucher MG07_7_10_clone2.*': 'An. albertoi',
        r'Anopheles strodei voucher MG07_7_10_clone3.*': 'An. albertoi',
        
        # An. strodei CPform
        r'Anopheles strodei s. l. CP form voucher ES20_4_1.*': 'An. strodei CPform',
        r'Anopheles strodei s. l. CP form voucher MG15_01_01.*': 'An. strodei CPform',
        r'Anopheles strodei s. l. CP form voucher MG15_06_12.*': 'An. strodei CPform',
        r'Anopheles strodei s. l. CP form voucher PR21_110.*': 'An. strodei CPform',
        
        # An. arthuri
        r'Anopheles arthuri voucher MG07_6_3.*': 'An. arthuri',
        r'Anopheles arthuri voucher FSP-USP:MG07_6_3.*': 'An. arthuri',
        r'Anopheles arthuri voucher MG24_1.*': 'An. arthuri',
        r'Anopheles arthuri voucher RO08_1.*': 'An. arthuri',
        r'Anopheles arthuri voucher RO08_104.*': 'An. arthuri',
        r'Anopheles arthuri voucher SP31_120.*': 'An. arthuri',
        
        # An. rondoni
        r'Anopheles rondoni voucher PR28_34_100.*': 'An. rondoni',
        r'Anopheles rondoni voucher PR28_36_02 clone 2.*': 'An. rondoni',
        r'Anopheles rondoni voucher PR28_36_02 clone 3.*': 'An. rondoni',
        r'Anopheles rondoni voucher PR28_36_02.*': 'An. rondoni',
        r'Anopheles rondoni voucher PR28_55_100 clone 2.*': 'An. rondoni',
        r'Anopheles rondoni voucher PR28_55_100.*': 'An. rondoni',
        
        # An. benarrochi
        r'Anopheles benarrochi voucher AC15_109.*': 'An. benarrochi',
        r'Anopheles benarrochi voucher AC18_115.*': 'An. benarrochi',
        r'Anopheles benarrochi voucher AC18_117.*': 'An. benarrochi',
        r'Anopheles benarrochi voucher AC18_120.*': 'An. benarrochi',
        
        # An. triannulatus
        r'Anopheles triannulatus voucher AC1_108.*': 'An. triannulatus',
        r'Anopheles triannulatus voucher AP17_04_01.*': 'An. triannulatus',
        r'Anopheles triannulatus voucher ES03_03_01.*': 'An. triannulatus',
        r'Anopheles triannulatus voucher MG56_12_03.*': 'An. triannulatus',
        r'Anopheles triannulatus voucher SP09_02.*': 'An. triannulatus',
        
        # An. oryzalimnetes
        r'Anopheles oryzalimnetes voucher SP09_03.*': 'An. oryzalimnetes',
        
        # An. deaneorum
        r'Anopheles deaneorum voucher AC01_07.*': 'An. deaneorum',
        r'Anopheles deaneorum voucher AC02_02.*': 'An. deaneorum',
        r'Anopheles deaneorum voucher MS08_127.*': 'An. deaneorum',
        
        # An. albitarsis
        r'Anopheles albitarsis voucher MG11_20_3.*': 'An. albitarsis',
        r'Anopheles albitarsis voucher SP104_2_2.*': 'An. albitarsis',
        r'Anopheles albitarsis voucher VP06_01_01.*': 'An. albitarsis',
        
        # An. marajoara
        r'Anopheles marajoara voucher AP21_50_1.*': 'An. marajoara',
        r'Anopheles marajoara voucher AP5_01_04.*': 'An. marajoara',
        r'Anopheles marajoara voucher PA3_1_13.*': 'An. marajoara',
        
        # An. braziliensis
        r'Anopheles braziliensis voucher AP21_39_3.*': 'An. braziliensis',
        r'Anopheles braziliensis voucher SP16_03.*': 'An. braziliensis',
        
        # An. argyritarsis
        r'Anopheles sawyeri voucher CE17_14_100 clone 1.*': 'An. argyritarsis',
        r'Anopheles sawyeri voucher CE17_14_100 clone 3.*': 'An. argyritarsis',
        r'Anopheles sawyeri voucher CE17_14_100 clone 5.*': 'An. argyritarsis',
        r'Anopheles sawyeri voucher CE17_14_100.*': 'An. argyritarsis',
        r'Anopheles argyritarsis voucher MG25_4.*': 'An. argyritarsis s.l.',
        
        # An. darlingi
        r'Anopheles darlingi voucher AC20_21_100.*': 'An. darlingi',
        r'Anopheles darlingi voucher AP13_03_06.*': 'An. darlingi',
        r'Anopheles darlingi voucher FSP-USP:AP13_03_06.*': 'An. darlingi',
        r'Anopheles darlingi voucher AP17_01_10.*': 'An. darlingi',
        r'Anopheles darlingi voucher FSP-USP:AP17_01_10.*': 'An. darlingi',
        
        # An. lanei
        r'Anopheles lanei voucher CJ02_02 clone 1.*': 'An. lanei',
        r'Anopheles lanei voucher CJ02_02 clone 2.*': 'An. lanei',
        r'Anopheles lanei voucher CJ02_02 clone 3.*': 'An. lanei',
        r'Anopheles lanei voucher CJ02_02.*': 'An. lanei',
        r'Anopheles lanei voucher CJ02_03.*': 'An. lanei',
        
        # An. atacamensis
        r'Anopheles atacamensis voucher FSP-USP:E-12992.*': 'An. atacamensis',
        r'Anopheles atacamensis voucher FSP-USP:E-12993.*': 'An. atacamensis',
        r'Anopheles atacamensis.*': 'An. atacamensis',
        
        # An. guarani
        r'Anopheles guarani voucher PR29 clone 1.*': 'An. guarani',
        r'Anopheles guarani voucher PR29.*': 'An. guarani',
        r'Anopheles guarani voucher PR29_08.*': 'An. guarani',
        r'Anopheles guarani voucher PR29_09_06.*': 'An. guarani',
        
        # An. parvus
        r'Anopheles parvus voucher AS5_1.*': 'An. parvus',
        r'Anopheles parvus voucher AS5_2.*': 'An. parvus',
        r'Anopheles parvus voucher AS5_3.*': 'An. parvus',
        r'Anopheles parvus voucher AS5_4.*': 'An. parvus',
        r'Anopheles parvus voucher MG07_9_1.*': 'An. parvus',
        r'Anopheles parvus voucher MG56_2.*': 'An. parvus',
        r'Anopheles parvus voucher MG07_9_20.*': 'An. parvus',
        r'Anopheles parvus voucher PR28_18_1.*': 'An. parvus',
        r'Anopheles parvus voucher PR28_5_1.*': 'An. parvus',
        r'Anopheles parvus voucher PR28_65_6.*': 'An. parvus',
        
        # An. antunesi
        r'Anopheles antunesi voucher RJ03_11.*': 'An. antunesi',
        r'Anopheles antunesi voucher RJ03_12.*': 'An. antunesi',
        r'Anopheles antunesi voucher RJ03_12_clone1.*': 'An. antunesi',
        r'Anopheles antunesi voucher RJ03_12_clone3.*': 'An. antunesi',
        r'Anopheles antunesi voucher RJ03_1_clone2.*': 'An. antunesi',
        r'Anopheles antunesi voucher RJ03_13.*': 'An. antunesi',
        r'Anopheles antunesi voucher RJ03_6.*': 'An. antunesi',
        r'Anopheles antunesi voucher VP09_17.*': 'An. antunesi',
        r'Anopheles antunesi voucher FSP-USP:VP11b.*': 'An. antunesi',
        r'Anopheles antunesi voucher VP11b.*': 'An. antunesi',
        
        # An. lutzii ss
        r'Anopheles lutzii voucher FSP-USP:SP02_9_2.*': 'An. lutzii ss',
        r'Anopheles lutzii voucher FSP-USP:SP02_10_-05.*': 'An. lutzii ss',
        r'Anopheles lutzii voucher SP02_11_9.*': 'An. lutzii ss',
        r'Anopheles lutzii voucher SP02_12_1.*': 'An. lutzii ss',
        r'Anopheles lutzii voucher SP02_13_3.*': 'An. lutzii ss',
        r'Anopheles lutzii voucher SP02_14_6.*': 'An. lutzii ss',
        r'Anopheles lutzii voucher SP02_15_5.*': 'An. lutzii ss',
        r'Anopheles lutzii voucher SP02_9_2.*': 'An. lutzii ss',
        r'Anopheles lutzii voucher SP02_10_05.*': 'An. lutzii ss',
        
        # An. lutzii B
        r'Anopheles lutzii voucher B369 clone 3.*': 'An. lutzii B',
        r'Anopheles lutzii voucher B369 clone 4.*': 'An. lutzii B',
        r'Anopheles lutzii voucher B369 clone 5.*': 'An. lutzii B',
        r'Anopheles lutzii s. l. 1 RS16 voucher RS16a.*': 'An. lutzii B',
        r'Anopheles lutzii s. l. 1 RS16 voucher RS16b.*': 'An. lutzii B',
        
        # An. lutzii A
        r'Anopheles lutzii voucher A325.*': 'An. lutzii A',
        r'Anopheles lutzii s. l. 2 RS19 voucher RS19_13.*': 'An. lutzii A',
        r'Anopheles lutzii s. l. 2 RS19 voucher RS19_21.*': 'An. lutzii A',
        r'Anopheles lutzii s. l. 2 RS19 voucher RS19_22.*': 'An. lutzii A',
        r'Anopheles lutzii s. l. 2 RS33 voucher RS33_105.*': 'An. lutzii A',
        r'Anopheles lutzii s. l. 2 RS33 voucher RS33a.*': 'An. lutzii A',
        r'Anopheles lutzii s. l. 2 RS33 voucher RS33b.*': 'An. lutzii A',
        
        # An. pristinus
        r'Anopheles pristinus voucher SP50a.*': 'An. pristinus',
        r'Anopheles pristinus voucher SP50b.*': 'An. pristinus',
        r'Anopheles pristinus voucher SP51_100.*': 'An. pristinus',
        r'Anopheles pristinus voucher SP53_100.*': 'An. pristinus',
        r'Anopheles pristinus voucher SP53_101.*': 'An. pristinus',
        r'Anopheles pristinus voucher SP53_4.*': 'An. pristinus',
        r'Anopheles pristinus voucher SP53_5.*': 'An. pristinus',
        r'Anopheles pristinus voucher SP55_2.*': 'An. pristinus',
        r'Anopheles pristinus voucher SP55_4.*': 'An. pristinus',
        r'Anopheles pristinus voucher VP11a.*': 'An. pristinus',
        
        # An. kompi
        r'Anopheles kompi voucher SP69_22_5.*': 'An. kompi',
        
        # An. cruzii
        r'Anopheles cruzii voucher ST16.*': 'An. cruzii'
    }
    
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



process_fasta('C:/Users/david/Downloads/filtered_coi.fasta', 'coi_agrupado.fasta')