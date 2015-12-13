if (($env:appveyor_repo_tag -eq "True") -and ($env:appveyor_repo_branch.StartsWith("v"))) {
    write-output "Deploying to PyPI..."
    Invoke-Expression "twine upload -u mcflugen -p $env:PYPI_PASS dist/*"
    write-output "OK"
} else {
    write-output "Not deploying."
    write-output "Repo tag is $env:appveyor_repo_tag"
    write-output "Branch name $env:appveyor_repo_branch"
}
