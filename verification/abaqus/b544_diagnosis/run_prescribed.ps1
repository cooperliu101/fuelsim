param(
    [string]$SourceDirectory = $PSScriptRoot,
    [string[]]$JobNames = @("b544_prescribed_default", "b544_prescribed_weak")
)

$ErrorActionPreference = "Stop"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
$ParentDirectory = Split-Path -Parent $SourceDirectory
$Work = Join-Path $env:TEMP ("fuelsim_b544_prescribed_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $ParentDirectory "extract_b544.py") $Work
foreach ($JobName in $JobNames) {
    Copy-Item (Join-Path $SourceDirectory "$JobName.inp") $Work
}

Push-Location $Work
try {
    foreach ($JobName in $JobNames) {
        & "C:\SIMULIA\Commands\abaqus.bat" job=$JobName input="$JobName.inp" cpus=1 output_precision=full interactive |
            Tee-Object -FilePath "${JobName}_run.log"
        if ($LASTEXITCODE -ne 0) { throw "Abaqus $JobName solve failed with exit code $LASTEXITCODE" }
        if (!(Test-Path "$JobName.sta") -or
            !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
            throw "Abaqus $JobName did not report successful completion"
        }
        $Outputs = @("${JobName}_nodal.csv", "${JobName}_integration.csv",
            "${JobName}_contact.csv", "${JobName}_energy.csv")
        & "C:\SIMULIA\Commands\abaqus.bat" python extract_b544.py "$JobName.odb" `
            $Outputs[0] $Outputs[1] $Outputs[2] $Outputs[3] |
            Tee-Object -FilePath "${JobName}_extract.log"
        if ($LASTEXITCODE -ne 0) { throw "Abaqus $JobName extraction failed with exit code $LASTEXITCODE" }
        if (!(Test-Path ($Outputs[0] + ".ok"))) { throw "Abaqus $JobName extraction did not complete" }
        foreach ($Output in $Outputs) {
            if (!(Test-Path $Output)) { throw "Abaqus extraction did not create $Output" }
            Copy-Item $Output $SourceDirectory
        }
        foreach ($Extension in @("dat", "msg", "sta")) {
            Copy-Item "$JobName.$Extension" $SourceDirectory
        }
        Copy-Item "${JobName}_run.log" $SourceDirectory
        if (Test-Path "${JobName}_extract.log") {
            Copy-Item "${JobName}_extract.log" $SourceDirectory
        }
        @("case=$JobName", "scope=prescribed native nodal U/T; no external loads; diagnostic only",
            "abaqus_version=3DEXPERIENCE R2018x", "cpus=1", "output_precision=full",
            "completed_utc=$([DateTime]::UtcNow.ToString('o'))", "work_directory=$Work",
            "input_sha256=$((Get-FileHash "$JobName.inp" -Algorithm SHA256).Hash)",
            "extractor_sha256=$((Get-FileHash 'extract_b544.py' -Algorithm SHA256).Hash)") +
            @($Outputs | ForEach-Object { "${_}_sha256=$((Get-FileHash $_ -Algorithm SHA256).Hash)" }) |
            Set-Content (Join-Path $SourceDirectory "${JobName}_provenance.txt")
        Write-Output "Completed $JobName in $Work"
    }
} finally {
    Pop-Location
}
