if (($env:appveyor_repo_tag -eq "True") -and ($env:appveyor_repo_branch.StartsWith("v"))) {
    Invoke-Expression "$env:PYTHON/Scripts/twine upload -u landlab -p $env:PYPI_PASS dist/*
}
