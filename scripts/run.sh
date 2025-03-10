prefix="antp"
data_path="../data"

echo 9 | gmx trjconv -f ${data_path}/md_OK_dt100.xtc \
                -o ${data_path}/pdb/${prefix}.pdb -sep \
                -dt 100

                # Il manque le .tpr