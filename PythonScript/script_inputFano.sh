for i in {20..168}
do

echo 1 25 25 ci_co_decay_sing ci_co_final_c_sing_$i ci_co_c_tricat_doub  fano_sing_c_doub$i\_
echo 1 25 25 ci_co_decay_sing ci_co_final_o_sing_$i ci_co_o_tricat_doub  fano_sing_o_doub$i\_

echo 1 35 25 ci_co_decay_trip ci_co_final_c_trip_$i ci_co_c_tricat_doub  fano_trip_c_doub$i\_
echo 1 35 25 ci_co_decay_trip ci_co_final_o_trip_$i ci_co_o_tricat_doub  fano_trip_o_doub$i\_

echo 1 35 10 ci_co_decay_trip ci_co_final_c_trip_$i ci_co_c_tricat_quad  fano_trip_c_quad$i\_
echo 1 35 10 ci_co_decay_trip ci_co_final_o_trip_$i ci_co_o_tricat_quad  fano_trip_o_quad$i\_

#echo 1 25 25 ci_co_decay ci_co_final_$i ci_co_dicat_doub fano$i\_

#echo 1 25 25 ci_co_decay ci_co_final_$i ci_co_dicat_doub fano$i\_

# echo 1 1 ci_anion_2 ci_anion_$i fano$i\_
# echo 1 25 10 ci_co_decay ci_co_final_$i ci_co_dicat_triplet  triplet_fano$i\_
# echo 1 9 6 ci_ch2_decay ci_ch2+_final_$i ci_ch2+_ion_singlet  singlet_fano$i\_
done
