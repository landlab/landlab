if (($env:appveyor_repo_tag -eq "true") -and ($env:appveyor_repo_tag_name.StartsWith("v"))) {
    write-output "Deploying to PyPI..."
    Invoke-Expression "twine upload -u mcflugen -p $env:PYPI_PASS dist/*"
    write-output "OK"
} else {
    write-output "Not deploying."
    write-output "Repo tag is $env:appveyor_repo_tag"
    write-output "Branch name $env:appveyor_repo_tag_name"
}
