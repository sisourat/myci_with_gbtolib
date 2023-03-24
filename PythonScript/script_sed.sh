declare -i imin=$1
declare -i imax=$2

ftemplate=$3
fout=$4
for i in $(seq $imin 1 $imax)
do
  echo $i
  sed "s/X/$i/g" $ftemplate > a
  python /home/nico/Workspace/GeneralFanoCICode/PythonScript/ormas_like.py a > b
  cat header.txt b > $fout$i
done
