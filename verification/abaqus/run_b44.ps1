param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_b44_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "extract_b44.py") $work
Set-Location $work

function Invoke-B44Case {
    param(
        [string]$Job,
        [string]$Mode,
        [string]$NodeFile,
        [string]$ContactFile,
        [int]$ExpectedNodeLines,
        [int]$ExpectedContactLines
    )
    $inputFile = $Job + ".inp"
    Copy-Item (Join-Path $SourceDirectory $inputFile) $work
    & "C:\SIMULIA\Commands\abaqus.bat" `
        "job=$Job" `
        "input=$inputFile" `
        output_precision=full `
        interactive
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus B4.4 $Mode solve failed with exit code $LASTEXITCODE"
    }
    if (!(Test-Path ($Job + ".sta")) -or
        !(Select-String -Path ($Job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus B4.4 $Mode solve did not report successful completion"
    }
    & "C:\SIMULIA\Commands\abaqus.bat" python extract_b44.py ($Job + ".odb") $NodeFile $ContactFile $Mode
    if ($LASTEXITCODE -ne 0 -or !(Test-Path $NodeFile) -or !(Test-Path $ContactFile) -or
        (Get-Content $NodeFile).Count -ne $ExpectedNodeLines -or
        (Get-Content $ContactFile).Count -ne $ExpectedContactLines) {
        throw "Abaqus B4.4 $Mode extraction failed"
    }
    Copy-Item $NodeFile $SourceDirectory
    Copy-Item $ContactFile $SourceDirectory
}

Invoke-B44Case "b44_hex8_nonmatching_partial_contact" "base" `
    "b44_hex8_nonmatching_partial_nodes.csv" "b44_hex8_nonmatching_partial_contact.csv" 589 145
Invoke-B44Case "b44_hex8_nonmatching_partial_contact_swapped" "swapped" `
    "b44_hex8_nonmatching_partial_swapped_nodes.csv" "b44_hex8_nonmatching_partial_swapped_contact.csv" 99 10
Write-Output "Abaqus B4.4 work directory: $work"
