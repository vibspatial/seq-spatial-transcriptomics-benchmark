# Exp66
pigz -dc -p 24 4789_Exp66_SC0009_S5_L001_R1_001.fastq.gz 4789_Exp66_SC0009_S5_L002_R1_001.fastq.gz | pigz -p 24 > 4789_Exp66_SC0009_S5_R1.fastq.gz
pigz -dc -p 24 4789_Exp66_SC0009_S5_L001_R2_001.fastq.gz 4789_Exp66_SC0009_S5_L002_R2_001.fastq.gz | pigz -p 24 > 4789_Exp66_SC0009_S5_R2.fastq.gz

# Exp74
pigz -dc -p 24 4789_Exp74_SC0010_S6_L001_R1_001.fastq.gz 4789_Exp74_SC0010_S6_L002_R1_001.fastq.gz | pigz -p 24 > 4789_Exp74_SC0010_S6_R1.fastq.gz
pigz -dc -p 24 4789_Exp74_SC0010_S6_L001_R2_001.fastq.gz 4789_Exp74_SC0010_S6_L002_R2_001.fastq.gz | pigz -p 24 > 4789_Exp74_SC0010_S6_R2.fastq.gz