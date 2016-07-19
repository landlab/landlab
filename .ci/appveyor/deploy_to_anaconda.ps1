if (($env:appveyor_repo_tag -eq "true") -and ($env:appveyor_repo_tag_name.StartsWith("v"))) {
    write-output "Deploying to anaconda main channel..."
    $channel = "main"
} else {
    write-output "Deploying to anaconda dev channel..."
    $channel = "dev"
    $env:BUILD_STR = "dev"
}

Invoke-Expression "conda config --set anaconda_upload no"

# write-output "Building package..."
# Invoke-Expression "conda build .conda -c landlab"

$file_to_upload = (conda build --output .conda) | Out-String
write-output "Uploading $file_to_upload..."
Invoke-Expression "anaconda -t $env:ANACONDA_TOKEN upload --force --user landlab --channel $channel $file_to_upload"
# Invoke-Expression "anaconda -t $env:ANACONDA_TOKEN upload --force --user landlab --channel $channel C:\\Miniconda\conda-bld\**\landlab*bz2"

write-output "OK"
