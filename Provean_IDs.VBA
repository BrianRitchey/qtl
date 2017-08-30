# Provean_IDs function in Excel VBA
#
# written by Brian Ritchey
#
#


Sub collect_PROVEAN_IDs()
‘ Queries PROVEAN website for AKR vs. DBA/2 aa changes (http://provean.jcvi.org/seq_submit.php).
‘ Collects PROVEAN urls into Excel cells for queries.
‘ When server is done processing, the data can be retrieved through these urls.

    Dim IE As New SHDocVw.InternetExplorer
    Dim HTMLDoc As MSHTML.HTMLDocument
    Dim HTMLInput_fasta As MSHTML.IHTMLElement
    Dim HTMLInput_variants As MSHTML.IHTMLElement
    Dim HTMLInput_inputtag As MSHTML.IHTMLElement
    Dim HTMLInput_inputtags As MSHTML.IHTMLElementCollection
    
    Dim genes As Range
    Dim gene As Range
    Dim query As String
    Dim var As Range
    Dim variants As Range
    Dim var_string As String
    Dim all_genes As Range
    Dim current_gene As Range
    Dim all_variants As Range
    Dim current_variant As Range
    Dim current_gene_offset As Range
    

    Set all_genes = Range("$A$1", Range("$A$1").End(xlDown))
    Set all_variants = Range("$D$1", Range("$D$1").End(xlDown))
   

    For Each current_gene In all_genes
    
        Set genes = Range(current_gene.Address, Range(current_gene.Address).End(xlToRight))
        
        Set current_gene_offset = Range(current_gene.Address)
        Set current_gene_offset = current_gene_offset.Offset(, 3)
        
        Set variants = Range(current_gene_offset.Address, Range(current_gene_offset.Address).End(xlToRight))
    For Each var In variants
      var_string = var_string & var.Cells & vbNewLine
    Next var
 query = genes.Item(2) & vbNewLine & genes.Item(3)
    
    IE.Visible = False
    IE.Navigate ("http://provean.jcvi.org/seq_submit.php")

    Do While IE.ReadyState <> READYSTATE_COMPLETE
    Loop
    
    Set HTMLDoc = IE.Document
    Set HTMLInput_fasta = HTMLDoc.getElementById("query")
    HTMLInput_fasta.Value = query
    Set HTMLInput_variants = HTMLDoc.getElementById("variant")
    HTMLInput_variants.Value = var_string
    Set HTMLInput_inputtags = HTMLDoc.getElementsByTagName("input")
    
'Submit the query
HTMLInput_inputtags(9).Click
query = ""
var_string = ""
Application.Wait (Now + TimeValue("0:00:10"))

'Create shell and get url which contains job id
For Each wd In CreateObject("Shell.Application").Windows
    If wd = "Internet Explorer" Then
        If InStr(wd.LocationName, "PROVEAN job submitted") <> 0 Then
            Exit For
        End If
    End If
Next wd

Dim Provean_ID As String
Provean_ID = wd.LocationURL

Dim row As Integer
row = genes.row

' change column number (15 in this case) based on output
Cells(row, 15).Value = Provean_ID
wd.Quit
Set wd = Nothing

Next current_gene

IE.Quit
    
End Sub

