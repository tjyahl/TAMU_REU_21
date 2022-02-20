import java.io.*;
import java.lang.*;


public class main {
        public static void main(String[] args) throws IOException{
		FileReader in = null;
		FileWriter out = null;
		try{
			in = new FileReader(args[0]);
			out = new FileWriter(args[1]);
		} catch(Exception e) {
			System.out.println(args[0]);
			System.out.println(args[1]);
			System.out.println("Wrong file");
			return;
		}

		
		int temp;
		String c = new String();
		boolean first = true;
		out.write('{');
		String graphString = new String();
		while((temp = in.read()) != -1){
			c = Character.toString(temp);
			if(c.equals("{")){
				graphString = new String();
				graphString += c;
				while(!(c = Character.toString(in.read())).equals("}")){
					graphString += c;
				}
				graphString += c;

				if(first){
					first = false;
				} else {
					out.write(',');
				}
				out.write(graphString.toCharArray());
			}
		}
		out.write('}');
		out.close();
        }
};

