package ch.irb.bamtools;

/**
 Copyright 2017 - Mathilde Foglierini Perez
 This code is distributed open source under the terms of the GNU Free Documention License.

 The Region object is a insert locus defined by its genomic coordinates and its gene Id (none if intergenic region)
 */
public class Region {
    private String chr;
    private int start;
    private int end;
    private String geneId;
    private String exonNumber;

    private boolean duplicate=false;

    public Region(String chr, int start, int end) {
        setChr(chr);
        setStart(start);
        setEnd(end);
    }

    public String getName() {
        return chr + ":" + start + "-" + end;
    }

    public int getStart() {
        return start;
    }

    public int getLength() {
        return (end - start + 1);
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public boolean isDuplicate() {
        return duplicate;
    }

    public void setDuplicate(boolean duplicate) {
        this.duplicate = duplicate;
    }
    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public String getGeneId() {
        return geneId;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public String getExonNumber() {
        return exonNumber;
    }

    public void setExonNumber(String exonNumber) {
        this.exonNumber = exonNumber;
    }
}
