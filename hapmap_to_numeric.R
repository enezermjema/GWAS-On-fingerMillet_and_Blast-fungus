##function to convert hapmap to numeric (-1,0 and 1)

## Assuming the initial hapmap data is a haploid formatted from TASSEL software

hapmap_numeric=function(hapmap_data="hap matrix"){
  
  ##pulls out first allele (ref) for a and second (alt) for b
  ref = substring(hapmap_data$alleles,1,1)
  
  #Same thing with the alt allele
  alt = substring(hapmap_data$alleles,3,3)
  
  #makes a copy of the hap matrix
  hap_num = hapmap_data
  
  #sets all allele values to NA
  hap_num[,12:ncol(hap_num)]=NA
  
  ## Turn allele a and allele b into 1 and -1.  Het into 0
  #line by line if a line is a then it places 1 in hap_num for the allele

  hap_num[x == a] = 1
  hap_num[x == b] = -1
  hap_num[x == "M"] = 0
  hap_num[x == "Y"] = 0
  hap_num[x == "K"] = 0
  hap_num[x == "R"] = 0
  hap_num[x == "W"] = 0
  hap_num[x == "S"] = 0
  hap_num[x== "N"]=NA
  
  return(hap_num)
}
```
