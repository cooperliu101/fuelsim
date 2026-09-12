param(
    [string[]]$Cases = @("gps_uniform_axial", "gps_uniform_axial_finite", "gps_two_slice_contact"),
    [string]$DestinationDirectory = $PSScriptRoot
)
$ErrorActionPreference = "Stop"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
New-Item -ItemType Directory -Force -Path $DestinationDirectory | Out-Null
foreach ($Case in ($Cases -join ',').Split(',')) {
    if ($Case -notin @("gps_uniform_axial", "gps_uniform_axial_finite", "gps_two_slice_contact")) {
        throw "Unexpected case: $Case"
    }
    $Work = Join-Path $env:TEMP ("fuelsim_gps_" + [guid]::NewGuid().ToString("N"))
    New-Item -ItemType Directory -Path $Work | Out-Null
    Copy-Item (Join-Path $PSScriptRoot "$Case.inp"), (Join-Path $PSScriptRoot "extract.py") $Work
    Write-Output "case=$Case work_directory=$Work"
    Push-Location $Work
    try {
        $ErrorActionPreference = "Continue"
        & C:\SIMULIA\Commands\abaqus.bat job=$Case input="$Case.inp" cpus=1 output_precision=full ask_delete=OFF interactive *> "${Case}_run.log"
        $AnalysisExitCode = $LASTEXITCODE
        $ErrorActionPreference = "Stop"
        foreach ($Extension in @("dat", "msg", "sta")) {
            if (Test-Path "$Case.$Extension") { Copy-Item "$Case.$Extension" $DestinationDirectory }
        }
        Get-Content "${Case}_run.log" | Set-Content (Join-Path $DestinationDirectory "${Case}_run.log") -Encoding UTF8
        if ($AnalysisExitCode -ne 0 -or !(Test-Path "$Case.sta") -or
            !(Select-String -Path "$Case.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
            throw "Abaqus failed; inspect $Work"
        }
        $ErrorActionPreference = "Continue"
        & C:\SIMULIA\Commands\abaqus.bat python extract.py "$Case.odb" $Case *> "${Case}_extract.log"
        $ExtractExitCode = $LASTEXITCODE
        $ErrorActionPreference = "Stop"
        Get-Content "${Case}_extract.log" | Set-Content (Join-Path $DestinationDirectory "${Case}_extract.log") -Encoding UTF8
        if ($ExtractExitCode -ne 0 -or !(Test-Path "${Case}_export_complete.txt") -or
            !(Select-String -Path "${Case}_export_complete.txt" -Pattern '^accepted_frames=10$' -Quiet)) {
            throw "Extraction failed or did not export every accepted increment; inspect $Work"
        }
        foreach ($Kind in @("nodes", "points", "contact")) {
            if (!(Test-Path "${Case}_${Kind}.csv")) { throw "Missing export ${Case}_${Kind}.csv" }
            Copy-Item "${Case}_${Kind}.csv" $DestinationDirectory
        }
        Copy-Item "${Case}_fields.txt" $DestinationDirectory
        @("case=$Case", "abaqus_version=3DEXPERIENCE R2018x", "cpus=1", "output_precision=full",
            "accepted_increments=10", "completed_utc=$([DateTime]::UtcNow.ToString('o'))", "work_directory=$Work",
            "input_sha256=$((Get-FileHash "$Case.inp" -Algorithm SHA256).Hash)",
            "extractor_sha256=$((Get-FileHash 'extract.py' -Algorithm SHA256).Hash)") |
            Set-Content (Join-Path $DestinationDirectory "${Case}_provenance.txt")
    } finally { Pop-Location }
}
