param([Parameter(Mandatory=$true)][string]$DestinationDirectory)
$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path $DestinationDirectory | Out-Null
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
foreach ($Case in @("expansion_small", "expansion_finite")) {
    $Work = Join-Path $env:TEMP ("fuelsim_expansion_" + [guid]::NewGuid().ToString("N"))
    New-Item -ItemType Directory -Path $Work | Out-Null
    Copy-Item (Join-Path $PSScriptRoot "$Case.inp") $Work
    Copy-Item (Join-Path $PSScriptRoot "../extract_b7_rz.py") $Work
    Push-Location $Work
    try {
        & C:\SIMULIA\Commands\abaqus.bat job=$Case input="$Case.inp" cpus=1 output_precision=full ask_delete=OFF interactive
        if (!(Test-Path "$Case.sta") -or !(Select-String -Path "$Case.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
            throw "Native calculation failed: $Work"
        }
        & C:\SIMULIA\Commands\abaqus.bat python extract_b7_rz.py "$Case.odb" $Case
        if ($LASTEXITCODE -ne 0) { throw "Extraction failed: $Work" }
        foreach ($Suffix in @(".sta", ".msg", ".dat", "_nodes.csv", "_points.csv")) {
            Copy-Item "$Case$Suffix" $DestinationDirectory
        }
        @("case=$Case", "abaqus_version=3DEXPERIENCE R2018x", "cpus=1", "output_precision=full",
          "completed_utc=$([DateTime]::UtcNow.ToString('o'))", "work_directory=$Work",
          "input_sha256=$((Get-FileHash "$Case.inp" -Algorithm SHA256).Hash)",
          "extractor_sha256=$((Get-FileHash extract_b7_rz.py -Algorithm SHA256).Hash)") |
          Set-Content (Join-Path $DestinationDirectory "${Case}_run.txt")
    } finally { Pop-Location }
}
