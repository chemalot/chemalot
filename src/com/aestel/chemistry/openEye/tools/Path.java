/*
   Copyright 2006-2014 Man-Ling Lee & Alberto Gobbi

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Contact: aestelSW@gmail.com
*/

package com.aestel.chemistry.openEye.tools;


/**
 * Non mutable path storing a sequence of nodes and vertices
 * @author albertgo
 *
 * @param <N>
 * @param <V>
 */
public class Path<N,V>
{  private V[] vertices;
   private N[] nodes;

   public Path(N[] nodes, V[] vertices)
   {  this.nodes = nodes;
      this.vertices = vertices;
   }
   
   /**
    * returns a safe copy of the node sequence.
    */
   public N[] getNodeSequence()
   {  return nodes.clone();
   }

   /**
    * returns a safe copy of the vertex sequence.
    */
   public V[] getVertexSequence()
   {  return vertices.clone();
   }
   
   public int getNNodes()     { return nodes.length; }
   public int getNVertices()  { return vertices.length; }
   public N   getLastNode()   { return nodes[nodes.length-1]; }
   public V   getLastVertix() { return vertices[vertices.length-1]; }
   
   
   public String toString()
   {  StringBuilder sb = new StringBuilder((vertices.length + nodes.length)*2);
      int nIdx = 0;
      for(N node : nodes)
      {  sb.append(node.toString());
         if(nIdx < vertices.length)
            sb.append(vertices[nIdx++].toString()); 
      }
      while(nIdx < vertices.length)
         sb.append(vertices[nIdx++].toString());
      
      return sb.toString();
   }
}
