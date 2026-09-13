$ErrorActionPreference="Stop"
$Source=$PSScriptRoot
& (Join-Path $Source "run_b61.ps1") -SourceDirectory $Source
& (Join-Path $Source "run_b526.ps1") -SourceDirectory $Source -Case "b526_friction_reversal"
