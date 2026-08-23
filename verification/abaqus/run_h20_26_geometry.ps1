param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_26_geometry_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
$jobs = @(
    @{ Name = "h20_26_hex20_sts_finite_sliding_probe"; Case = "finite_sliding" },
    @{ Name = "h20_26_hex20_sts_tilted_normal_probe"; Case = "tilted_normal" }
)
Copy-Item (Join-Path $SourceDirectory "extract_h20_26_geometry.py") $work
foreach ($job in $jobs) {
    Copy-Item (Join-Path $SourceDirectory ($job.Name + ".inp")) $work
}
Set-Location $work

foreach ($job in $jobs) {
    & "C:\SIMULIA\Commands\abaqus.bat" `
        job=$($job.Name) `
        input=$($job.Name + ".inp") `
        output_precision=full `
        interactive
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus H20.26 geometry probe failed with exit code $LASTEXITCODE"
    }
    if (!(Test-Path ($job.Name + ".sta")) -or
        !(Select-String -Path ($job.Name + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus H20.26 geometry probe did not report successful completion: $($job.Name)"
    }
    $nodal = $job.Name + "_nodal.csv"
    $history = $job.Name + "_history.csv"
    & "C:\SIMULIA\Commands\abaqus.bat" `
        python extract_h20_26_geometry.py `
        ($job.Name + ".odb") `
        $job.Case `
        $nodal `
        $history
    if (!(Test-Path $nodal) -or !(Test-Path $history)) {
        throw "Abaqus H20.26 geometry extraction did not create both CSV files: $($job.Name)"
    }
    Copy-Item $nodal $SourceDirectory
    Copy-Item $history $SourceDirectory
}
Write-Output "Abaqus H20.26 geometry-probe work directory: $work"
