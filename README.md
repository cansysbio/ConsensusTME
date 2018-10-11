# Consensus<sup>TME</sup> - Tumour microenvironment cell estimation 


`ConsensusTME` is a integrative tool for R that uses an consensus approach to generating cancer specific signatures for multiple cell types found within the tumour microenvironment. 

These consensus gene sets are then used within a ssGSEA framework to provide normalised enrichment scores for each of the cell types representing the relative abundance of cell types across multiple samples. 

#
<p align="center">
  <img src="https://github.com/cansysbio/ConsensusTME/blob/master/Overview.png" width="520" height="412"></div>
</p>

# Installation

A package for Consensus<sup>TME</sup> is currently in development 

# Usage

-Until an R package is available consensus gene sets have been made available that can be used in various statistical frameworks including ssGSEA. 

-It is important to note that results generated in this way will provide quantification of cell types that are relative across samples rather than across cell types.

-Consensus<sup>TME</sup> gene sets have been made to be cancer specific, if analysis is aimed at healthy tissue or PBMCs better results may be obtained by using the unspecific gene sets

-Consensus gene sets can be generated manually through cloning of respository N.B. Paths within functions will need to be updated.
