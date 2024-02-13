# **Genome Size Estimation with Jellyfish**
* memory decided via hash size (size of the genome plus seq error) T. dalmanni genome is 0.62 gb, so i'll assume 0.7 (nanoplot doesnt accept less than 1G so do that)
* not sure what kmer size to use, 21 seems sensible 
* jellyfish V 2.2.10

Ran [jellyfish.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/jellyfish.sh) on each species

### **kmer distribution plots**
* used [genomescope](http://qb.cshl.edu/genomescope/) to visualise the .histo files generated by jellyfish (read length=10000, max kmer coverage=1000000):
- [_T. whitei_ Female](http://genomescope.org/analysis.php?code=blR2SdZ6dlrwedM2Fgs8)
- [_T. whitei_ Male]
<br><br>
* All seem acceptable, containing two peaks
* However, est. genome size seems odd, being variable, and larger in females than males in _D. meigenii_ and _T. whitei_ and the opposite in _T. dalmanni_