'''
Created on May 26, 2017

@author: eze
'''
from io import StringIO
import unittest

from SNDG.Comparative.VcfSnpeffIO import VcfSnpeffIO


class TestVcfSnpeffIO(unittest.TestCase):

    
    missense = u"""
    NC_000962.3     14785   .       T       C       2020.0  .       AC=1;AF=1.00;AN=1;DP=77;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=60.00;QD=26.23;SOR=0.772;ANN=C|missense_variant|MODERATE|Rv0012|Rv0012|transcript|Rv0012|protein_coding|1/1|c.697T>C|p.Cys233Arg|697/789|697/789|233/262||,C|upstream_gene_variant|MODIFIER|Rv0008c|Rv0008c|transcript|Rv0008c|protein_coding||c.-2474A>G|||||2474|,C|upstream_gene_variant|MODIFIER|Rv0010c|Rv0010c|transcript|Rv0010c|protein_coding||c.-1227A>G|||||1227|,C|upstream_gene_variant|MODIFIER|Rv0011c|Rv0011c|transcript|Rv0011c|protein_coding||c.-790A>G|||||790|,C|upstream_gene_variant|MODIFIER|trpG|Rv0013|transcript|Rv0013|protein_coding||c.-129T>C|||||129|,C|downstream_gene_variant|MODIFIER|gyrA|Rv0006|transcript|Rv0006|protein_coding||c.*4967T>C|||||4967|,C|downstream_gene_variant|MODIFIER|Rv0007|Rv0007|transcript|Rv0007|protein_coding||c.*3957T>C|||||3957|WARNING_TRANSCRIPT_NO_START_CODON,C|downstream_gene_variant|MODIFIER|ppiA|Rv0009|transcript|Rv0009|protein_coding||c.*1769T>C|||||1769|,C|downstream_gene_variant|MODIFIER|pknB|Rv0014c|transcript|Rv0014c|protein_coding||c.*805A>G|||||805|,C|downstream_gene_variant|MODIFIER|pknA|Rv0015c|transcript|Rv0015c|protein_coding||c.*2682A>G|||||2682|,C|downstream_gene_variant|MODIFIER|pbpA|Rv0016c|transcript|Rv0016c|protein_coding||c.*3974A>G|||||3974| GT:AD:DP:GQ:PL  1:0,77:77:99:2050,0
    """
    
    rna = u"""
    NC_000962.3    1472362    .    C    T    2123.0    .    AC=1;AF=1.00;AN=1;BaseQRankSum=-1.588;ClippingRankSum=0.000;DP=92;FS=6.426;MLEAC=1;MLEAF=1.00;MQ=59.92;MQRankSum=3.830;QD=23.08;ReadPosRankSum=1.277;SOR=0.148;ANN=T|upstream_gene_variant|MODIFIER|Rv1313c|Rv1313c|transcript|Rv1313c|protein_coding||c.-2857G>A|||||2857|WARNING_TRANSCRIPT_NO_START_CODON,T|upstream_gene_variant|MODIFIER|Rv1314c|Rv1314c|transcript|Rv1314c|protein_coding||c.-2110G>A|||||2110|,T|downstream_gene_variant|MODIFIER|atpC|Rv1311|transcript|Rv1311|protein_coding||c.*4682C>T|||||4682|,T|downstream_gene_variant|MODIFIER|Rv1312|Rv1312|transcript|Rv1312|protein_coding||c.*4231C>T|||||4231|,T|downstream_gene_variant|MODIFIER|murA|Rv1315|transcript|Rv1315|protein_coding||c.*785C>T|||||785|WARNING_TRANSCRIPT_NO_START_CODON,T|downstream_gene_variant|MODIFIER|ogt|Rv1316c|transcript|Rv1316c|protein_coding||c.*4772G>A|||||4772|,T|intragenic_variant|MODIFIER|rrs|Rvnr01|gene_variant|Rvnr01|||n.1472362C>T||||||    GT:AD:DP:GQ:PL    1:6,86:92:99:2153,0
    """
    
    upstream = u"""
    NC_000962.3     34044   .       T       C       1657.0  .       AC=1;AF=1.00;AN=1;DP=57;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=60.00;QD=29.59;SOR=0.841;ANN=C|upstream_gene_variant|MODIFIER|bioF2|Rv0032|transcript|Rv0032|protein_coding||c.-251T>C|||||251|,C|upstream_gene_variant|MODIFIER|acpA|Rv0033|transcript|Rv0033|protein_coding||c.-2563T>C|||||2563|,C|upstream_gene_variant|MODIFIER|Rv0034|Rv0034|transcript|Rv0034|protein_coding||c.-2823T>C|||||2823|,C|upstream_gene_variant|MODIFIER|fadD34|Rv0035|transcript|Rv0035|protein_coding||c.-3215T>C|||||3215|,C|downstream_gene_variant|MODIFIER|Rv0024|Rv0024|transcript|Rv0024|protein_coding||c.*4837T>C|||||4837|WARNING_TRANSCRIPT_NO_START_CODON,C|downstream_gene_variant|MODIFIER|Rv0025|Rv0025|transcript|Rv0025|protein_coding||c.*4437T>C|||||4437|,C|downstream_gene_variant|MODIFIER|Rv0026|Rv0026|transcript|Rv0026|protein_coding||c.*2976T>C|||||2976|WARNING_TRANSCRIPT_NO_START_CODON,C|downstream_gene_variant|MODIFIER|Rv0027|Rv0027|transcript|Rv0027|protein_coding||c.*2538T>C|||||2538|,C|downstream_gene_variant|MODIFIER|Rv0028|Rv0028|transcript|Rv0028|protein_coding||c.*2225T>C|||||2225|,C|downstream_gene_variant|MODIFIER|Rv0029|Rv0029|transcript|Rv0029|protein_coding||c.*890T>C|||||890|WARNING_TRANSCRIPT_NO_START_CODON,C|downstream_gene_variant|MODIFIER|Rv0030|Rv0030|transcript|Rv0030|protein_coding||c.*491T>C|||||491|,C|intergenic_region|MODIFIER|Rv0030-bioF2|Rv0030-Rv0032|intergenic_region|Rv0030-Rv0032|||n.34044T>C||||||        GT:AD:DP:GQ:PL  1:0,56:56:99:1687,0
    """
    
    frameshift = u"""
    NC_000962.3    4359165    .    G    GTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC    1264.97    .    AC=1;AF=1.00;AN=1;DP=27;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=68.11;QD=26.86;SOR=0.693;ANN=GTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC|frameshift_variant|HIGH|espK|Rv3879c|transcript|Rv3879c|protein_coding|1/1|c.617_618insGCCGGGAACCCCGATCACCCCGGGAACCCCGATCA|p.Pro219fs|617/2190|617/2190|206/729||,GTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC|downstream_gene_variant|MODIFIER|espI|Rv3876|transcript|Rv3876|protein_coding||c.*4155_*4156insTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC|||||4156|,GTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC|downstream_gene_variant|MODIFIER|eccD1|Rv3877|transcript|Rv3877|protein_coding||c.*2623_*2624insTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC|||||2624|,GTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC|downstream_gene_variant|MODIFIER|espJ|Rv3878|transcript|Rv3878|protein_coding||c.*1630_*1631insTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC|||||1631|,GTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC|downstream_gene_variant|MODIFIER|espL|Rv3880c|transcript|Rv3880c|protein_coding||c.*1033_*1034insGCCGGGAACCCCGATCACCCCGGGAACCCCGATCA|||||1033|WARNING_TRANSCRIPT_NO_START_CODON,GTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC|downstream_gene_variant|MODIFIER|espB|Rv3881c|transcript|Rv3881c|protein_coding||c.*1377_*1378insGCCGGGAACCCCGATCACCCCGGGAACCCCGATCA|||||1377|,GTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC|downstream_gene_variant|MODIFIER|eccE1|Rv3882c|transcript|Rv3882c|protein_coding||c.*2866_*2867insGCCGGGAACCCCGATCACCCCGGGAACCCCGATCA|||||2866|,GTGATCGGGGTTCCCGGGGTGATCGGGGTTCCCGGC|downstream_gene_variant|MODIFIER|mycP1|Rv3883c|transcript|Rv3883c|protein_coding||c.*4251_*4252insGCCGGGAACCCCGATCACCCCGGGAACCCCGATCA|||||4251|WARNING_TRANSCRIPT_NO_START_CODON;LOF=(espK|Rv3879c|1|1.00)    GT:AD:DP:GQ:PL    1:0,23:23:99:1304,0
    """
    
    stop_gained = u"""
    NC_000962.3    4280467    .    C    A    2006.0    .    AC=1;AF=1.00;AN=1;BaseQRankSum=-1.647;ClippingRankSum=0.000;DP=75;FS=3.069;MLEAC=1;MLEAF=1.00;MQ=60.00;MQRankSum=0.000;QD=26.75;ReadPosRankSum=-0.375;SOR=0.308;ANN=A|stop_gained|HIGH|Rv3815c|Rv3815c|transcript|Rv3815c|protein_coding|1/1|c.322G>T|p.Glu108*|322/756|322/756|108/251||,A|upstream_gene_variant|MODIFIER|Rv3813c|Rv3813c|transcript|Rv3813c|protein_coding||c.-1252G>T|||||1252|,A|upstream_gene_variant|MODIFIER|Rv3814c|Rv3814c|transcript|Rv3814c|protein_coding||c.-452G>T|||||452|,A|upstream_gene_variant|MODIFIER|Rv3817|Rv3817|transcript|Rv3817|protein_coding||c.-1180C>A|||||1180|WARNING_TRANSCRIPT_NO_START_CODON,A|upstream_gene_variant|MODIFIER|Rv3818|Rv3818|transcript|Rv3818|protein_coding||c.-1982C>A|||||1982|WARNING_TRANSCRIPT_NO_START_CODON,A|upstream_gene_variant|MODIFIER|Rv3819|Rv3819|transcript|Rv3819|protein_coding||c.-3529C>A|||||3529|,A|downstream_gene_variant|MODIFIER|csp|Rv3811|transcript|Rv3811|protein_coding||c.*4050C>A|||||4050|,A|downstream_gene_variant|MODIFIER|PE_PGRS62|Rv3812|transcript|Rv3812|protein_coding||c.*2382C>A|||||2382|WARNING_TRANSCRIPT_NO_START_CODON,A|downstream_gene_variant|MODIFIER|Rv3816c|Rv3816c|transcript|Rv3816c|protein_coding||c.*325G>T|||||325|WARNING_TRANSCRIPT_NO_START_CODON,A|downstream_gene_variant|MODIFIER|papA2|Rv3820c|transcript|Rv3820c|protein_coding||c.*3952G>T|||||3952|WARNING_TRANSCRIPT_NO_START_CODON    GT:AD:DP:GQ:PL    1:1,74:75:99:2036,0
    """
    
    insertion = u"""
    NC_000962.3    3932774    .    C    CGGCGCCGGCGGCGCGGGTGGCGCCGGCGCGGACGCGACCGCTACCGGCGCCACCGGCGGCACCGGGTTCGCCGGT    78.97    .    AC=1;AF=1.00;AN=1;DP=5;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=57.47;QD=26.32;SOR=2.833;ANN=CGGCGCCGGCGGCGCGGGTGGCGCCGGCGCGGACGCGACCGCTACCGGCGCCACCGGCGGCACCGGGTTCGCCGGT|disruptive_inframe_insertion|MODERATE|PE_PGRS54|Rv3508|transcript|Rv3508|protein_coding|1/1|c.1781_1782insCGCGGGTGGCGCCGGCGCGGACGCGACCGCTACCGGCGCCACCGGCGGCACCGGGTTCGCCGGTGGCGCCGGCGG|p.Gly594_Ala595insAlaGlyGlyAlaGlyAlaAspAlaThrAlaThrGlyAlaThrGlyGlyThrGlyPheAlaGlyGlyAlaGlyGly|1782/5706|1782/5706|594/1901||INFO_REALIGN_3_PRIME,CGGCGCCGGCGGCGCGGGTGGCGCCGGCGCGGACGCGACCGCTACCGGCGCCACCGGCGGCACCGGGTTCGCCGGT|downstream_gene_variant|MODIFIER|PE_PGRS53|Rv3507|transcript|Rv3507|protein_coding||c.*2060_*2061insGGCGCCGGCGGCGCGGGTGGCGCCGGCGCGGACGCGACCGCTACCGGCGCCACCGGCGGCACCGGGTTCGCCGGT|||||2061|,CGGCGCCGGCGGCGCGGGTGGCGCCGGCGCGGACGCGACCGCTACCGGCGCCACCGGCGGCACCGGGTTCGCCGGT|downstream_gene_variant|MODIFIER|ilvX|Rv3509c|transcript|Rv3509c|protein_coding||c.*4102_*4103insACCGGCGAACCCGGTGCCGCCGGTGGCGCCGGTAGCGGTCGCGTCCGCGCCGGCGCCACCCGCGCCGCCGGCGCC|||||4102|WARNING_TRANSCRIPT_NO_START_CODON    GT:AD:DP:GQ:PL    1:0,3:3:99:118,0
    """

    deletion = u"""
    NC_000962.3    4189134    .    TGCGCCTACA    T    1086.97    .    AC=1;AF=1.00;AN=1;DP=25;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=60.00;QD=29.02;SOR=1.609;ANN=T|disruptive_inframe_deletion|MODERATE|Rv3737|Rv3737|transcript|Rv3737|protein_coding|1/1|c.1439_1447delGCCTACAGC|p.Arg480_Gln482del|1439/1590|1439/1590|480/529||INFO_REALIGN_3_PRIME,T|upstream_gene_variant|MODIFIER|Rv3733c|Rv3733c|transcript|Rv3733c|protein_coding||c.-4631_-4623delTGTAGGCGC|||||4623|,T|upstream_gene_variant|MODIFIER|tgs2|Rv3734c|transcript|Rv3734c|protein_coding||c.-3253_-3245delTGTAGGCGC|||||3245|,T|downstream_gene_variant|MODIFIER|Rv3735|Rv3735|transcript|Rv3735|protein_coding||c.*2558_*2566delGCGCCTACA|||||2558|,T|downstream_gene_variant|MODIFIER|Rv3736|Rv3736|transcript|Rv3736|protein_coding||c.*1440_*1448delGCGCCTACA|||||1440|,T|downstream_gene_variant|MODIFIER|PPE66|Rv3738c|transcript|Rv3738c|protein_coding||c.*142_*150delTGTAGGCGC|||||150|,T|downstream_gene_variant|MODIFIER|PPE67|Rv3739c|transcript|Rv3739c|protein_coding||c.*1141_*1149delTGTAGGCGC|||||1149|,T|downstream_gene_variant|MODIFIER|Rv3740c|Rv3740c|transcript|Rv3740c|protein_coding||c.*1690_*1698delTGTAGGCGC|||||1698|,T|downstream_gene_variant|MODIFIER|Rv3741c|Rv3741c|transcript|Rv3741c|protein_coding||c.*3036_*3044delTGTAGGCGC|||||3044|,T|downstream_gene_variant|MODIFIER|Rv3742c|Rv3742c|transcript|Rv3742c|protein_coding||c.*3707_*3715delTGTAGGCGC|||||3715|WARNING_TRANSCRIPT_NO_START_CODON,T|downstream_gene_variant|MODIFIER|ctpJ|Rv3743c|transcript|Rv3743c|protein_coding||c.*4248_*4256delTGTAGGCGC|||||4256|WARNING_TRANSCRIPT_NO_START_CODON    GT:AD:DP:GQ:PL    1:0,25:25:99:1126,0
    """
    
    not_annotated = u"""
    NC_000962.3     4198611 .       CG      C       2689.97 .       AC=1;AF=1.00;AN=1;DP=84;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=60.13;QD=32.02;SOR=0.843 GT:AD:DP:GQ:PL  1:0,84:84:99:2729,0
    """
    
    VCF_HEADER = u"""
            #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  1109_S1_L00
            """

    

    def test_parse_all(self):
        for x in [TestVcfSnpeffIO.missense,TestVcfSnpeffIO.upstream,
                  TestVcfSnpeffIO.frameshift,TestVcfSnpeffIO.stop_gained,
                  TestVcfSnpeffIO.insertion,TestVcfSnpeffIO.deletion,
                  TestVcfSnpeffIO.not_annotated]:
             
            variants = list(VcfSnpeffIO.parse( StringIO(TestVcfSnpeffIO. VCF_HEADER +  x) ))
            assert variants
     
    def test_missense(self):
        variant,effects = list(VcfSnpeffIO.parse( 
            StringIO(TestVcfSnpeffIO. VCF_HEADER +  TestVcfSnpeffIO.missense) ))[0]
        effect = effects[0]
        self.assertEqual(14785,variant.POS)
        self.assertEqual(233,effect.aa_pos)
        self.assertEqual("C",effect.aa_ref)
        self.assertEqual("R",effect.aa_alt)
    
    def test_rna_var(self):
        _,effects = list(VcfSnpeffIO.parse( 
            StringIO(TestVcfSnpeffIO. VCF_HEADER +  TestVcfSnpeffIO.rna) ))[0]
        self.assertEqual("intragenic_variant",effects[0].annotation[0])
        
    
    def test_insertion(self):
        _,effects = list(VcfSnpeffIO.parse( 
            StringIO(TestVcfSnpeffIO. VCF_HEADER +  TestVcfSnpeffIO.insertion) ))[0]
        effect = effects[0]
        
        #Gly594_Ala595insAlaGlyGlyAlaGlyAlaAspAlaThrAlaThrGlyAlaThrGlyGlyThrGlyPheAlaGlyGlyAlaGlyGly
        
        self.assertEqual(594,effect.aa_pos)
        self.assertEqual("G",effect.aa_ref) 
            # it only asks for G and not GA because if the hgvs_p is an interval, snpeff shows an interval and not the original aa Seq
        self.assertEqual("AGGAGADATATGATGGTGFAGGAGG",effect.aa_alt)
    
    def test_deletion(self):
        _,effects = list(VcfSnpeffIO.parse( 
            StringIO(TestVcfSnpeffIO. VCF_HEADER +  TestVcfSnpeffIO.deletion) ))[0]
        effect = effects[0]        
        #Arg480_Gln482del        
        self.assertEqual(480,effect.aa_pos)
        self.assertEqual("R",effect.aa_ref) 
            # it only asks for R and not RE because if the hgvs_p is an interval, snpeff shows an interval and not the original aa Seq
        self.assertEqual("del",effect.aa_alt)
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestVcfSnpeffIO.testName']
    unittest.main()