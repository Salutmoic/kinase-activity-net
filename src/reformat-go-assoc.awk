function join(array, start, end, sep, result, i)
{
    if (sep == "")
       sep = " "
    else if (sep == SUBSEP) # magic value
       sep = ""
    result = array[start]
    for (i = start + 1; i <= end; i++)
        result = result sep array[i]
    return result
}

{
    # Only take experimentally verified annotations
    if ($3 != "EXP" && $3 != "IDA" && $3 != "IPI" &&
            $3 != "IMP" && $3 != "IGI" && $3 != "IEP")
        next
    if(a[$1])
    {
        a[$1] = a[$1]";"$2;
    }
    else
    {
        a[$1]=$2;
    }
}
END {
    OFS = "\t"
    for (i in a)
    {
        if (index(a[i], ";"))
        {
            split(a[i], terms, ";")
            asort(terms)
            termlist = join(terms, 1, length(terms), ";")
        }
        else
        {
            termlist = a[i] 
        }
        print i, termlist;
    }
}
