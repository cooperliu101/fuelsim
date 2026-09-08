param([Parameter(Mandatory = $true)][string]$SourceDirectory,
      [Parameter(Mandatory = $true)][string]$WorkDirectory)
$ErrorActionPreference = 'Stop'
$JobName = 'c3d20rt_finite_contact_frictionless'
$SourceDirectory = (Resolve-Path $SourceDirectory).Path
$WorkDirectory = (Resolve-Path $WorkDirectory).Path
Write-Output "reference_work_directory=$WorkDirectory"
if ((Get-FileHash (Join-Path $SourceDirectory "$JobName.inp")).Hash -ne
    (Get-FileHash (Join-Path $WorkDirectory "$JobName.inp")).Hash) {
    throw 'Native contact frames require the matching reference input'
}
Copy-Item (Join-Path $SourceDirectory 'extract_c3d20rt_contact_frames.py') $WorkDirectory
Set-Location $WorkDirectory
& 'C:\SIMULIA\Commands\abaqus.bat' python extract_c3d20rt_contact_frames.py "$JobName.odb" "${JobName}_frames.csv"
if ($LASTEXITCODE -ne 0) { throw 'Native contact frame extraction failed' }
if ((Get-Content "${JobName}_frames.csv").Count -ne 3641) {
    throw 'Expected 280 complete contact frames with 13 nodes each'
}
Copy-Item "${JobName}_frames.csv" $SourceDirectory
Write-Output 'contact_frame_rows=3640'
Write-Output 'extraction_completed=true'
