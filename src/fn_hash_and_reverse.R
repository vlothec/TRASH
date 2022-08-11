

Hash_And_Reverse = function(repeat_seq, kmer = 4)
{
  #for leading strand
  {
    repeat_test = repeat_seq
    #divide sequence into kmers
    #hash kmers
    repeat_test = strsplit(repeat_test, split = "")[[1]]
    temp_repeat_test = repeat_test
    for(i in 1 : length(repeat_test))
    {
      if(i <= (length(repeat_test) - kmer + 1))
      {
        repeat_test[i] = paste(temp_repeat_test[i:(i+kmer-1)], collapse = "")
      } else
      {
        repeat_test[i] = paste(paste(temp_repeat_test[i:length(temp_repeat_test)], collapse = ""), paste(temp_repeat_test[1:(kmer-(length(temp_repeat_test)-i)-1)], collapse = ""), sep = "")
      }
    }
    hash_val = vector(mode = "numeric", length = length(repeat_test))
    for(i in 1 : length(repeat_test))
    {
      hash_single = strsplit(repeat_test[i], split = "")[[1]]
      for(j in 1 : kmer)
      {
        actg_val = 0 * (hash_single[j] == "A") + 1 * (hash_single[j] == "C") + 2 * (hash_single[j] == "T") + 3 * (hash_single[j] == "G") 
        hash_val[i] = hash_val[i] + (4 ^ j) * actg_val
      }
    }
    #for all possible start positions
    lowest_hash_start = sum(hash_val * 1:length(hash_val))
    lowest_hash_id = 1
    for(i in 2 : length(repeat_test))
    {
      #multiply hash values by kmer position and extract the sum value
      hash_val = c(hash_val[2: length(repeat_test)], hash_val[1])
      hash_sum = sum(hash_val * 1:length(hash_val))
      if(hash_sum < lowest_hash_start)
      {
        lowest_hash_id = i
        lowest_hash_start = hash_sum
      }
    }
    lowest_hash_leading = lowest_hash_start
    lowest_hash_id_leading = lowest_hash_id
  }
  remove(lowest_hash_id, lowest_hash_start, hash_val, repeat_test, temp_repeat_test)
  #for rev comp
  {
    {
      repeat_test = revCompString(repeat_seq)
      #divide sequence into kmers
      #hash kmers
      repeat_test = strsplit(repeat_test, split = "")[[1]]
      temp_repeat_test = repeat_test
      for(i in 1 : length(repeat_test))
      {
        if(i <= (length(repeat_test) - kmer + 1))
        {
          repeat_test[i] = paste(temp_repeat_test[i:(i+kmer-1)], collapse = "")
        } else
        {
          repeat_test[i] = paste(paste(temp_repeat_test[i:length(temp_repeat_test)], collapse = ""), paste(temp_repeat_test[1:(kmer-(length(temp_repeat_test)-i)-1)], collapse = ""), sep = "")
        }
      }
      hash_val = vector(mode = "numeric", length = length(repeat_test))
      for(i in 1 : length(repeat_test))
      {
        hash_single = strsplit(repeat_test[i], split = "")[[1]]
        for(j in 1 : kmer)
        {
          actg_val = 0 * (hash_single[j] == "A") + 1 * (hash_single[j] == "C") + 2 * (hash_single[j] == "T") + 3 * (hash_single[j] == "G") 
          hash_val[i] = hash_val[i] + (4 ^ j) * actg_val
        }
      }
      #for all possible start positions
      lowest_hash_start = sum(hash_val * 1:length(hash_val))
      lowest_hash_id = 1
      for(i in 2 : length(repeat_test))
      {
        #multiply hash values by kmer position and extract the sum value
        hash_val = c(hash_val[2: length(repeat_test)], hash_val[1])
        hash_sum = sum(hash_val * 1:length(hash_val))
        if(hash_sum < lowest_hash_start)
        {
          lowest_hash_id = i
          lowest_hash_start = hash_sum
        }
      }
      lowest_hash_rev_comp = lowest_hash_start
      lowest_hash_id_rev_comp = lowest_hash_id
    }
  }
  remove(lowest_hash_id, lowest_hash_start, hash_val, repeat_test, temp_repeat_test)
  
  #find lowest hash between leading and rev comp
  if(lowest_hash_rev_comp <= lowest_hash_leading)
  {
    temp_repeat_test = strsplit(repeat_seq, split = "")[[1]]
    #find lowest hash sum and change the repeat to that orientation
    if(lowest_hash_id_leading == 1)
    {
      repeat_fixed = temp_repeat_test
    } else
    {
      repeat_fixed = c(temp_repeat_test[lowest_hash_id_leading : length(temp_repeat_test)], temp_repeat_test[1 : (lowest_hash_id_leading - 1)])
    }
  } else
  {
    temp_repeat_test = strsplit(revCompString(repeat_seq), split = "")[[1]]
    #find lowest hash sum and change the repeat to that orientation
    if(lowest_hash_id_rev_comp == 1)
    {
      repeat_fixed = temp_repeat_test
    } else
    {
      repeat_fixed = c(temp_repeat_test[lowest_hash_id_rev_comp : length(temp_repeat_test)], temp_repeat_test[1 : (lowest_hash_id_rev_comp - 1)])
    }
  }
  
  repeat_fixed = paste(repeat_fixed, collapse = "")
  return(repeat_fixed)
}

