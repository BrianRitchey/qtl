# Navigate_to_PROVEAN Excel VBA subroutine
#
# written by Brian Ritchey
#
#


Sub Navigate_to_PROVEAN()
' Need to navigate to PROVEAN url in order to get results.
‘ This subroutine navigates to PROVEAN query webpages for all queries.

Dim IE As New SHDocVw.InternetExplorer
Dim HTMLDoc As MSHTML.HTMLDocument

Dim provean_urls As Range
'change this to appropriate cell
Set provean_urls = Range("O1", Range("O1").End(xlDown))
Dim url As Range
Dim http As String
    IE.Visible = False
For Each url In provean_urls
    http = url
    IE.Navigate (http)
    Do While IE.ReadyState <> READYSTATE_COMPLETE
    Loop
Next url

End Sub
