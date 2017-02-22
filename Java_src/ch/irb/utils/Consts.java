package ch.irb.utils;

import java.io.File;
import java.util.*;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 Copyright 2017 - Mathilde Foglierini Perez
 This code is distributed open source under the terms of the GNU Free Documention License.

Utils class.
 */

public final class Consts {

	public static final String ls = System.getProperty("line.separator");
	public static final String fs = System.getProperty("file.separator");
	public static final String[] aminoAcidsList = { "A", "G", "K", "N", "S", "D", "V", "L", "T", "I", "E", "Q",
			"R", "H", "M", "P", "C", "F", "W", "Y", "X", "-" };
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public static final ArrayList<String> aminoAcids = new ArrayList(Arrays.asList(aminoAcidsList));
	
	@SuppressWarnings("rawtypes")
	public static <K extends Comparable, V extends Comparable> Map<K, V> sortByValuesAsc(Map<K, V> map) {
		List<Map.Entry<K, V>> entries = new LinkedList<Map.Entry<K, V>>(map.entrySet());

		Collections.sort(entries, new Comparator<Map.Entry<K, V>>() {

			@SuppressWarnings("unchecked")
			public int compare(Entry<K, V> o1, Entry<K, V> o2) {
				return o1.getValue().compareTo(o2.getValue());
			}
		});

		// LinkedHashMap will keep the keys in the order they are inserted
		// which is currently sorted on natural ordering
		Map<K, V> sortedMap = new LinkedHashMap<K, V>();

		for (Map.Entry<K, V> entry : entries) {
			sortedMap.put(entry.getKey(), entry.getValue());
		}

		return sortedMap;
	}

	@SuppressWarnings("rawtypes")
	public static <K extends Comparable, V extends Comparable> Map<K, V> sortByValuesDesc(Map<K, V> map) {
		List<Map.Entry<K, V>> entries = new LinkedList<Map.Entry<K, V>>(map.entrySet());

		Collections.sort(entries, new Comparator<Map.Entry<K, V>>() {

			@SuppressWarnings("unchecked")
			public int compare(Entry<K, V> o1, Entry<K, V> o2) {
				return o2.getValue().compareTo(o1.getValue());
			}
		});

		// LinkedHashMap will keep the keys in the order they are inserted
		// which is currently sorted on natural ordering
		Map<K, V> sortedMap = new LinkedHashMap<K, V>();

		for (Map.Entry<K, V> entry : entries) {
			sortedMap.put(entry.getKey(), entry.getValue());
		}

		return sortedMap;
	}
	
    /*
     * Get the extension of a file.
     */  
    public static String getExtension(File f) {
        String ext = null;
        String s = f.getName();
        int i = s.lastIndexOf('.');

        if (i > 0 &&  i < s.length() - 1) {
            ext = s.substring(i+1).toLowerCase();
        }
        return ext;
    }	

}
