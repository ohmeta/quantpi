a="--output a.gz"
array=(${a// / })
echo ${array[@]}[1]

