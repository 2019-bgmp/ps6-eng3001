#!/bin/python
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=nt_counter   ### Job Name
#SBATCH --output=nt.out         ### File in which to store job output
#SBATCH --error=Fnt.err         ### File in which to store job error messages
#SBATCH --time=0-00:60:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp          ### Account used for job submission

FILE1="/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_1"
FILE2="/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_2"
FILE3="/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq.unmatched"

NT_sum1 = 0
NT_sum2 = 0
NT_sum3 = 0
avg_len_array = []

#Loop through all 3 files to count the number of NTs in each file
with open(FILE1, "r") as fh:
    LN = 0 #Line number
    for line in fh:
        LN+=1
        line = line.strip('\n') #Manipulate the sequence line
        if LN % 4 == 2:
            length = len(line)
            NT_sum1+=length
            avg_len_array.append(length)

with open(FILE2, "r") as fh:
    LN = 0 #Line number
    for line in fh:
        LN+=1
        line = line.strip('\n') #Manipulate the sequence line
        if LN % 4 == 2:
            length=0
            length = len(line)
            NT_sum2+=length
            avg_len_array.append(length)

with open(FILE3, "r") as fh:
    LN = 0 #Line number
    for line in fh:
        LN+=1
        line = line.strip('\n') #Manipulate the sequence line
        if LN % 4 == 2:
            length=0
            length = len(line)
            NT_sum3+=length
            avg_len_array.append(length)

#calculating the total sum of all the reads
Total_sum = NT_sum1 + NT_sum2 + NT_sum3
#calculating the coverage
Coverage = Total_sum / 2000000.0
#calculation of average length of every read in all three files
avg_len = float(sum(avg_len_array))/float(len(avg_len_array))
print ("num records:", len(avg_len_array))


print("File1 length: ", NT_sum1)
print("File2 length: ", NT_sum2)
print("File3 length: ", NT_sum3)
print("Total length: ", Total_sum)
print("Coverage: ", float(Coverage))
print("Average Length:", avg_len)

#calculating k-mer coverage for 31,41, and 49
k_mer_cov31 = (Coverage * (avg_len - 31 + 1)) / avg_len
k_mer_cov41 = (Coverage * (avg_len - 41 + 1)) / avg_len
k_mer_cov49 = (Coverage * (avg_len - 49 + 1)) / avg_len

print("K-mer coverage 31:", k_mer_cov31)
print("K-mer coverage 41:", k_mer_cov41)
print("K-mer coverage 49:", k_mer_cov49)
