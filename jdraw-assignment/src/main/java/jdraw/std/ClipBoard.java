package jdraw.std;
import java.util.LinkedList;
import java.util.List;

import jdraw.framework.*;

public final class ClipBoard {
	private static List<Figure> figures = new LinkedList<Figure>();
	
	public ClipBoard(){
	}
	
	public static List<Figure> get(){
		return new LinkedList<Figure>(figures);
	}
	
	public static void add(Figure f){
		figures.add(f);
	}
	
	public static void clear(){
		figures.clear();
	}
	
	public static void remove(Figure f){
		figures.remove(f);
	}
}
