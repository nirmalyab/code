import java.util.*;
import java.io.*;
import java.util.regex.*;

public class RelProg {
	
	BufferedReader relFile= null;
	// Create hashmap for left and right neighbors
	HashMap<Integer, HashSet<Integer>> leftNeigh = null;
	HashMap<Integer, HashSet<Integer>> rightNeigh = null;
	
	public RelProg(String relFileName)
	{
		try {
			relFile = new BufferedReader(new FileReader(relFileName));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		leftNeigh = new HashMap<Integer, HashSet<Integer>>();
		rightNeigh = new HashMap<Integer, HashSet<Integer>>();
		
		String regex = "^(\\S+)\\s+(\\S+)";
		Pattern pat = Pattern.compile(regex);
		Matcher mat = null;
		
		String line = null;
		
		try {
			while(null != (line = relFile.readLine()))
			{
				mat = pat.matcher(line);
				
				if(true == mat.find())
				{
					int left = Integer.parseInt(mat.group(1));
					int right = Integer.parseInt(mat.group(2));
					
					if(false ==  leftNeigh.containsKey(right))
					{
						HashSet<Integer> temp = new HashSet<Integer>();
						leftNeigh.put(right, temp);
					}
					
					HashSet<Integer> temp2 = leftNeigh.get(right);
					
					temp2.add(left);
					
					// Similarly do it for the rightNeigh
					
					if(false == rightNeigh.containsKey(left))
					{
						HashSet<Integer> temp = new HashSet<Integer>();
						rightNeigh.put(left, temp);
					
					}
					HashSet<Integer> temp3 = rightNeigh.get(left);
					temp3.add(right);
					
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			relFile.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
		
		
		
		
	public int [] getLeftNeigh(int gene)
	{
		Set<Integer> neigh = leftNeigh.get(gene);
		int [] neighArr = null;
		
		if(null != neigh)
		{
		
		neighArr = new int[neigh.size()];
		}
		else
		{
			neighArr = new int[0];
		}
		
		if(0 < neighArr.length)
		{
			Iterator<Integer> it = neigh.iterator();
			
			int i =0;
			while(true == it.hasNext())
			{
				int nextVar = it.next();
				
				neighArr[i] = nextVar;
				
				i++;
			}
		}
		
		return neighArr;
	}
	
	
	public int getLeftNeighCount(int gene)
	{
		Set<Integer> neigh = leftNeigh.get(gene);
		if(null == neigh)
			return 0;
		else
			return neigh.size();
	}
	
	public int [] getRightNeigh(int gene)
	{
		Set<Integer> neigh = rightNeigh.get(gene);
		
		int [] neighArr = null;
		if(null != neigh)
		{
		
		neighArr = new int[neigh.size()];
		}
		else
		{
			neighArr = new int[0];
		}
		
		if(0 < neighArr.length)
		{
			Iterator<Integer> it = neigh.iterator();
			
			int i =0;
			while(true == it.hasNext())
			{
				int nextVar = it.next();
				
				neighArr[i] = nextVar;
				
				i++;
			}
		}
		
		return neighArr;
	}
	
	public int getRightNeighCount(int gene)
	{
		Set<Integer> neigh = rightNeigh.get(gene);
		if(null == neigh)
			return 0;
		else
			return neigh.size();
	}
	
	public static void main(String []args)
	{
		RelProg rp = new RelProg(args[0]);
	}
	
}


