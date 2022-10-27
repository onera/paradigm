#
# Exemple of use :
#  - ./script/prepare_ci/run_top_test_base.sh --log_file_name tot
#  - ./script/prepare_ci/run_top_test_base.sh --log_file_name tot -k pdm_t_part_to_block_geom
#

while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        echo $1 $2 # // Optional to see the parameter:value result
   fi
   shift
done


execute_test () {
  # J'affiche son nom et demande confirmation pour l'effacer
  test_cmd=$*
  eval $*
  status=$?
  if [ "${status}" -ne "0" ]; then
    echo "${test_cmd} : Failed" >> ${log_file_name}
    echo "Failed"
    # test -f paradigm_0.log && echo pdm_t_nuga_adapt_cells && exit
    if [ "$(find -name '*core*' | wc -l)" ]; then
      rm *core*
    fi
  else
    echo "${test_cmd} : OK" >> ${log_file_name}
    echo "OK"
  fi
}


if [ -n "$log_file_name" ]; then
  echo "Log in file " $log_file_name
else
  log_file_name="output_test_log"
fi
rm $log_file_name

# echo $k
if [ -n "$k" ]; then
  echo "Execute test : " $k
  source script/prepare_ci/run_tests_$k.sh
else
  list=("pdm_t_dcube_nodal_gen" "pdm_t_part_to_block_geom")
  for item in ${list[*]}; do
    echo $item;
    source script/prepare_ci/run_tests_$item.sh
  done
fi
