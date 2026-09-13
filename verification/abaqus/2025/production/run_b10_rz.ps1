param(
    [Parameter(Mandatory = $true)][string[]]$Cases,
    [Parameter(Mandatory = $true)][string]$DestinationDirectory
)
function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path $DestinationDirectory | Out-Null
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
foreach ($Case in ($Cases -join ',').Split(',')) {
    if ($Case -notmatch '^b1[0-9]+_cax8t_[a-z_]+$') { throw "Unexpected case: $Case" }
    $Work = Join-Path $env:TEMP ("fuelsim_cax8t_" + [guid]::NewGuid().ToString("N"))
    New-Item -ItemType Directory -Path $Work | Out-Null
    Copy-Item (Join-Path $PSScriptRoot "$Case.inp") $Work
    Copy-Item (Join-Path $PSScriptRoot "b10*mesh.inc") $Work -ErrorAction SilentlyContinue
    $Extractor = "extract_b10_rz.py"
    if ($Case -match "contact") { $Extractor = "extract_b102_rz.py" }
    if ($Case -match "friction|sliding|recovery") { $Extractor = "extract_b8_rz.py" }
    Copy-Item (Join-Path $PSScriptRoot $Extractor) $Work
    Write-Output "case=$Case work_directory=$Work"
    Push-Location $Work
    try {
        Invoke-Abaqus job=$Case input="$Case.inp" output_precision=full ask_delete=OFF cpus=1 interactive
        foreach ($Extension in @("dat", "msg", "sta")) {
            if (Test-Path "$Case.$Extension") { Copy-Item "$Case.$Extension" $DestinationDirectory }
        }
        if (!(Test-Path "$Case.sta") -or
            !(Select-String -Path "$Case.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
            throw "Abaqus failed; inspect $Work"
        }
        Invoke-Abaqus python $Extractor "$Case.odb" $Case
        if ($LASTEXITCODE -ne 0) { throw "Extraction failed: $Case" }
        $Kinds = @("nodes", "points")
        if ($Case -match "friction|sliding|recovery") { $Kinds += "contact" }
        foreach ($Kind in $Kinds) {
            $Output = "${Case}_${Kind}.csv"
            if (!(Test-Path $Output) -or (Get-Content $Output | Measure-Object -Line).Lines -lt 2) {
                throw "Missing extracted field: $Output"
            }
            Copy-Item $Output $DestinationDirectory
        }
        @("case=$Case", "abaqus_version=Abaqus 2025 RELr427", "", "output_precision=full",
            "completed_utc=$([DateTime]::UtcNow.ToString('o'))", "work_directory=$Work",
            "input_sha256=$((Get-FileHash "$Case.inp" -Algorithm SHA256).Hash)",
            "extractor_sha256=$((Get-FileHash $Extractor -Algorithm SHA256).Hash)") |
            Set-Content (Join-Path $DestinationDirectory "${Case}_run.txt")
    } finally { Pop-Location }
}
