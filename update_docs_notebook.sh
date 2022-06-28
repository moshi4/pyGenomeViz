# Update docs notebook with lastest pyGenomeViz
for notebook_file in `find ./docs -name '*.ipynb' -type f`
do
    jupyter-nbconvert \
    --ClearMetadataPreprocessor.enabled=True \
    --execute --inplace --to notebook $notebook_file
done
