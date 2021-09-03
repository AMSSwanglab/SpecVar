sort -k1,1 $1 -o $1
cat $1 | awk -F"\t" 'BEGIN{
       pre_key=""
}{
       key=$1
       if(NR==1){
               pre_key=$1
               field2=$2
       }else if(key!="" && key==pre_key){
               pre_key=$1
               field2=field2","$2                           

       }else if (key!="" && key!=pre_key){
               print pre_key"\t"field2
               pre_key=$1
               field2=$2
       }

}END{
               print pre_key"\t"field2
       }' > ${1}_output
\mv ${1}_output ${1}
