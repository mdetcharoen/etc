1. ทดลองผลของการใช้ปุ๋ยสามชนิดต่อความยาวของใบหญ้าคากลุ่มละ 16 ใบได้ผลความยาว (cm) ของใบดังตาราง จงสร้างกราฟจาก R ggplot2 เพื่อนำเสนอข้อมูลนี้				
				
| leaf#  | 	control | 	fertilizer1 | 	fertilizer2	 | fertilizer3
| :---: | :---: | :---: | :---: | :---: |
1	 | 4.5 | 	6	 | 8.5 | 	4.5
2	 | 	4.9	 | 	6.8	 | 	9	 | 	5.5
3		 | 5.8	 | 	6.7	 | 	6	 | 6
4	 | 	6.1		 | 3.4	 | 	6.8	 | 	7.8
5		 | 6.4	 | 	6.9	 | 	7.8	 | 	6.3
6	 | 	8	 | 	7.1	 | 	7.6	 | 	6.2
7		 | 5.2	 | 	7.5	 | 	7.1	 | 	5.4
8		 | 4.5	 | 	6.2	 | 	6.1	 | 	4.4
9		 | 4.3	 | 	4.8	 | 	5	 | 	5.2
10	 | 	3.8	 | 	5.9	 | 	9.8	 | 	5.6
11	 | 	4.6	 | 	6.5	 | 	9.2	 | 	4.2
12	 | 	7.5	 | 	6.8	 | 	10.8	 | 	4.1
13	 | 	6.2	 | 	6.9	 | 	6.8	 | 	4.8
14	 | 	6.1	 | 	7	 | 	7.5	 | 	6.8
15	 | 	5.4	 | 	7.8	 | 	7.6	 | 	9.8
16	 | 	5.1	 | 	7.9		 | 7.8	 | 	7.1
</br>

				
2. เมื่อใช้คำสั่งในโปรแกรม fastp เป็น `fastp -i SRR26633039_1.fastq -o SRR26633039_1_qc.fastq -I SRR26633039_2.fastq -O SRR26633039_2_qc.fastq` อยากทราบว่าข้อมูลนี้เป็นข้อมูล single end หรือ paired end
</br>


4. ปกติแล้วค่าคุณภาพของเบสที่ได้จากการ sequencing จะถูกแสดงในรูปแบบของค่า Q เช่น Q20 และ Q30 อยากทราบว่าไฟล์ fastq ที่มีค่า Q20 = 94% นั้นค่า Q20 คืออะไร	
</br>


5. ในการทำ <i>de novo</i> genome assembly ของแบคทีเรียชนิดหนึ่ง เมื่อได้ assembly graph อยู่ในรูปไฟล์ `.gfa` แล้วนำไปเปิดดูพบว่าไม่ได้มีลักษณะเป็นวงตามที่คาดหมาย เมื่อลองเปลี่ยน parameters ต่าง ๆ ใน genome assembler แล้วก็ไม่ได้เป็นวง จงอิบายสาเหตุที่เป้นไปได้และข้อเสนอแนะในการปรับปรุง

6. จากข้อมูลคำสั่งที่กำหนด นักศึกษาตอบว่าผู้ใช้กำลังทำงานใดพร้อมบอกรายละเอียดของคำสั่งที่ใช้
   ```
   fasterq-dump SRR26633039
   fastp -h
   fastp -i SRR26633039_1.fastq -o SRR26633039_1_qc.fastq -I SRR26633039_2.fastq -O SRR26633039_2_qc.fastq
   spades.py
   spades.py --isolate -1 SRR26633039_1_qc.fastq -2 SRR26633039_2_qc.fastq --only-assembler -o outspades
   fasterq-dump SRR32385919
   fastp -i SRR32385919.fastq -o SRR32385919_qc.fastq
   flye -h
   flye --nano-hq SRR32385919_qc.fastq -t 14 -o flyeout
   minimap2 -ax map-ont -t 14 flyeout/selected_sequences.fasta SRR32385919_qc.fastq > mapped.sam
   samtools view -b mapped.sam > mapped.bam
   samtools sort mapped.bam > mapped.sorted.bam
   samtools index mapped.sorted.bam
   pilon --help
   pilon --genome flyeout/selected_sequences.fasta --bam mapped.sorted.bam --output pilon_ --changes
 ```

