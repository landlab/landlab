if (($env:appveyor_repo_tag -eq "True") -and ($env:appveyor_repo_branch.StartsWith("v"))) {
    Invoke-Expression "twine upload -u mcflugen -p $env:PYPI_PASS dist/*"
}
