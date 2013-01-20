     BEGIN {FS=";"}                                             
     {                                                          
       if (NR % 2 == 1) offset=$13;                             
       else print $0 "loadtime=;" offset;                      
     }                                                          
