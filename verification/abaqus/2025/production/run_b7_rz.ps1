param(
    [Parameter(Mandatory = $true)][string[]]$Cases,
    [Parameter(Mandatory = $true)][string]$DestinationDirectory,
    [string]$SourceDirectory
)
function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
if ([string]::IsNullOrWhiteSpace($SourceDirectory)) { $SourceDirectory = $PSScriptRoot }
New-Item -ItemType Directory -Force -Path $DestinationDirectory | Out-Null
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
foreach ($Case in ($Cases -join ',').Split(',')) {
    if ($Case -notmatch '^b7[0-9]+_rz_[a-z_]+$') { throw "Unexpected case name $Case" }
    $Work = Join-Path $env:TEMP ("fuelsim_rz_" + [guid]::NewGuid().ToString("N"))
    New-Item -ItemType Directory -Path $Work | Out-Null
    Copy-Item (Join-Path $SourceDirectory "$Case.inp") $Work
    Copy-Item (Join-Path $SourceDirectory "b7_rz_*mesh.inc") $Work
    Copy-Item (Join-Path $SourceDirectory "extract_b7_rz.py") $Work
    Write-Output "case=$Case work_directory=$Work"
    $Arguments = '/d /c start "" /b /wait "' + $env:ComSpec +
        '" /d /c C:\SIMULIA\Commands\abq2025.bat' +
        " job=$Case input=$Case.inp output_precision=full ask_delete=OFF cpus=1 interactive"
    $Watch = [System.Diagnostics.Stopwatch]::StartNew()
    $Process = Start-Process $env:ComSpec -ArgumentList $Arguments -WorkingDirectory $Work -PassThru -Wait
    $Watch.Stop()
    foreach ($Extension in @("dat", "msg", "sta")) {
        if (Test-Path (Join-Path $Work "$Case.$Extension")) {
            Copy-Item (Join-Path $Work "$Case.$Extension") $DestinationDirectory
        }
    }
    if ($Process.ExitCode -ne 0 -or !(Test-Path (Join-Path $Work "$Case.sta")) -or
        !(Select-String -Path (Join-Path $Work "$Case.sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus failed: $Case; inspect $Work"
    }
    Push-Location $Work
    try {
        Invoke-Abaqus python extract_b7_rz.py "$Case.odb" $Case
        if ($LASTEXITCODE -ne 0) { throw "Extraction failed: $Case" }
        foreach ($Output in @("${Case}_nodes.csv", "${Case}_points.csv")) {
            if (!(Test-Path $Output) -or (Get-Content $Output | Measure-Object -Line).Lines -lt 2) {
                throw "Extraction produced no data: $Output"
            }
        }
        Copy-Item "${Case}_nodes.csv", "${Case}_points.csv" $DestinationDirectory
    } finally { Pop-Location }
    @("case=$Case", "external_seconds=$($Watch.Elapsed.TotalSeconds)", "", "cpus=1",
      "input_sha256=$((Get-FileHash (Join-Path $Work "$Case.inp") -Algorithm SHA256).Hash)") |
        Set-Content (Join-Path $DestinationDirectory "${Case}_run.txt")
}
