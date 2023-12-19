
Rscript readextract4QC_module2.R

if [ $? -eq 0 ]; then
echo "runs well"
download_success=true
else
  echo "run failed"
download_success=false
fi
# 检查下载状态
while [ "$download_success" = false ]; do
Rscript readextract4QC_module2.R
if [ $? -eq 0 ]; then
echo "runs well"
download_success=true
else
  echo "run failed"
download_success=false
fi
done
