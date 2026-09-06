param(
    [Parameter(Mandatory = $true)][string[]]$Cases,
    [Parameter(Mandatory = $true)][string]$DestinationDirectory
)
$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path $DestinationDirectory | Out-Null
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
foreach ($Case in ($Cases -join ',').Split(',')) {
    if ($Case -notmatch '^b8[0-9]_rz_[a-z_]+$') { throw "Unexpected case: $Case" }
    $Work = Join-Path $env:TEMP ("fuelsim_rz_sliding_" + [guid]::NewGuid().ToString("N"))
    New-Item -ItemType Directory -Path $Work | Out-Null
    Copy-Item (Join-Path $PSScriptRoot "$Case.inp"), (Join-Path $PSScriptRoot "b8_rz_sliding_mesh.inc"),
        (Join-Path $PSScriptRoot "extract_b8_rz.py") $Work
    Write-Output "case=$Case work_directory=$Work"
    Push-Location $Work
    try {
        & C:\SIMULIA\Commands\abaqus.bat job=$Case input="$Case.inp" cpus=1 output_precision=full ask_delete=OFF interactive
        foreach ($Extension in @("dat", "msg", "sta")) {
            if (Test-Path "$Case.$Extension") { Copy-Item "$Case.$Extension" $DestinationDirectory }
        }
        if (!(Test-Path "$Case.sta") -or
            !(Select-String -Path "$Case.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
            throw "Abaqus failed; inspect $Work"
        }
        & C:\SIMULIA\Commands\abaqus.bat python extract_b8_rz.py "$Case.odb" $Case
        if ($LASTEXITCODE -ne 0) { throw "Extraction failed: $Case" }
        foreach ($Kind in @("nodes", "points", "contact")) {
            $Output = "${Case}_${Kind}.csv"
            if (!(Test-Path $Output) -or (Get-Content $Output | Measure-Object -Line).Lines -lt 2) {
                throw "Missing extracted field: $Output"
            }
            Copy-Item $Output $DestinationDirectory
        }
        @("case=$Case", "abaqus_version=3DEXPERIENCE R2018x", "cpus=1", "output_precision=full",
            "completed_utc=$([DateTime]::UtcNow.ToString('o'))", "work_directory=$Work",
            "input_sha256=$((Get-FileHash "$Case.inp" -Algorithm SHA256).Hash)",
            "mesh_sha256=$((Get-FileHash 'b8_rz_sliding_mesh.inc' -Algorithm SHA256).Hash)",
            "extractor_sha256=$((Get-FileHash 'extract_b8_rz.py' -Algorithm SHA256).Hash)") |
            Set-Content (Join-Path $DestinationDirectory "${Case}_run.txt")
    } finally { Pop-Location }
}
